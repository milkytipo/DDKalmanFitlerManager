
/**
 * Created by wuzida on 2019/5/13.
 * Notification: there maybe is existing a problem that  the conversion between primitive type(double type operation)
 * that it may generate a possible loss of precision. Float and double type are used to reduce storage space when
 storing rather than keeping high accuracy, frequent conversion data type operation may lead related the result to a big error.
 */
public class DDKalmanFilterManager {
    private double Ac[][] = {{1,0},{0,1}};
    private static int Sp = 1;
    private static double Sv=0.5;
    private static double Pxyz0[][] ={{30,0},{0,3}};
    private static double[] Rv =  {0.5,0.1};
    private double Qp0[][] = {{Sp+Sv/3, Sv/2},{Sv/2,Sv}};
    private double nxtState_pos[];
    private double nxtState_vel[];
    private double Pk[][] = new double[6][6];
    private double Qw[][] = new double[6][6];
    private double  satpos_actv[][];
    private double  satvel_actv[][] ;
    private double  satpos_ref_actv[][];
    private double  satvel_ref_actv[][] ;
    private double  satpos_rot_corr[][]; //storing the sat positions after earth rotation corrections
    private double obs_actv[];
    private  double obs_ref_actv[];
    private double obs_corr_actv[] ;
    private double[] deltaP_DD;
    private double[] sat_vel_DD;
    private double[] bur = new double[3];

    private double transmitTime_actv[] ;
    private double satClkCorr_actv[][] ;
    private double satClkCorr_ref_actv[][] ;
    private double cn0_actv[] ;
    private double  az[] ;
    private double  el[] ;
    private double deltaP[];
    private double deltaP_ref[];
    private double trop[];
    private double iono[] ;
    private double DOP[] ;
    private double H[][] ;
    private double H_DD[][] ;
    public  static double[] pos_ref_xyz= {-2853445.926, 4667466.476, 3268291.272};
    public static double   C = 299792458;
    public static double PI =3.141592653589793;
    //TODO:maybe such definitaion would cause unexpected trouble
    public  static class  PvtCalculator{
        public double[] stt_x  = new double[2];  //stt_* refer the bur,namely baseline's length
        public double[] stt_y = new double[2];
        public double[] stt_z  = new double[2];
        public double P[][] = new double[6][6];  //convirance
    };

    public DDKalmanFilterManager(double[] x,double[] y, double[] z) //xyz contain the position and velocity of user
    {
        PvtCalculator pvtCalculator  = new PvtCalculator();
        pvtCalculator.stt_x[0] = x[0]-pos_ref_xyz[0];
        pvtCalculator.stt_y[0] = y[0]-pos_ref_xyz[1];
        pvtCalculator.stt_z[0] = z[0]-pos_ref_xyz[2];
        pvtCalculator.stt_x[1] = x[0];
        pvtCalculator.stt_y[1] = y[1];
        pvtCalculator.stt_z[1] = z[2];
        pvtCalculator.P = MatrixBlock(Pxyz0,3);
        this.Qw =MatrixBlock(Qp0,3);
    }
    public DDKalmanFilterManager() //JUST FOR TEST
    {
        PvtCalculator pvtCalculator  = new PvtCalculator();

        pvtCalculator.P = MatrixBlock(Pxyz0,3);
        this.Qw =MatrixBlock(Qp0,3);
    }

    public void Update(int[] activeChannel,double[][] satpos,double[][] satpos_ref,double[] obs,double[] obs_ref, PvtCalculator pvtCalculator,  double[] doppSmooth, double[] doppSmooth_ref, double[] elfore,double[] azfore,double[]cn0)
    {
        int  nmbOfSatellites = activeChannel.length;
        satpos_actv =  new double[3][nmbOfSatellites];
        satvel_actv=  new double[3][nmbOfSatellites];
        satpos_ref_actv =  new double[3][nmbOfSatellites];
        satvel_ref_actv=  new double[3][nmbOfSatellites];
        satpos_rot_corr =  new double[3][nmbOfSatellites]; //storing the sat positions after earth rotation corrections
        deltaP = new double[nmbOfSatellites];
        deltaP_ref = new double[nmbOfSatellites];
        obs_actv = new double[nmbOfSatellites];
        obs_ref_actv = new double[nmbOfSatellites];
        cn0_actv= new double[nmbOfSatellites];
        az = new double[nmbOfSatellites];
        el = new double[nmbOfSatellites];;
        //       double DOP= new double[5];

        for(int i =0;i<nmbOfSatellites;i++)
        {
            // System.arraycopy(satpos_ref[i], 0, aa[i], 0, a[0].length);
            for(int j = 0;j<3;j++)
            {
                satpos_actv[j][i] = satpos[j][activeChannel[i]-1];
                satvel_actv[j][i]= satpos[j+3][activeChannel[i]-1];
                satpos_ref_actv[j][i] = satpos_ref[j][activeChannel[i]-1];
                satvel_ref_actv[j][i]= satpos_ref[j+3][activeChannel[i]-1];
            }
            deltaP[i]   = doppSmooth[activeChannel[i]-1];
            deltaP_ref[i]   = doppSmooth_ref[activeChannel[i]-1];
            obs_actv[i] = obs[activeChannel[i]-1];
            obs_ref_actv[i]=obs_ref[activeChannel[i]-1];
            cn0_actv[i] = cn0[ activeChannel[i]-1];
            el[i]   = elfore[activeChannel[i]-1];
            az[i]   = azfore[activeChannel[i]-1];
        }
        int max_el = find(el,1);
        PvtCalculator KalFilt =  pvtCalculator;
        nxtState_pos = new double[]{KalFilt.stt_x[0], KalFilt.stt_y[0], KalFilt.stt_z[0]};
        nxtState_vel = new double[]{KalFilt.stt_x[1], KalFilt.stt_y[1], KalFilt.stt_z[1]};
        for(int i =0;i<bur.length;i++){
            bur[i] =nxtState_pos[i] ;
        }
        for (int i = 0;i<3;i++){
            nxtState_pos[i] = nxtState_pos[i] + pos_ref_xyz[i];
        }
        double[] LLH = new double[3];
        LLH = ECEF2geo( nxtState_pos[0], nxtState_pos[1], nxtState_pos[2] );

        double [][] sat2usr_mtrx = new double[3][nmbOfSatellites];
        for(int i =0;i<nmbOfSatellites;i++) {
            sat2usr_mtrx[0][i]  =   satpos_rot_corr[0][i] - pos_ref_xyz[0];
            sat2usr_mtrx[1][i]  =   satpos_rot_corr[1][i] - pos_ref_xyz[1];
            sat2usr_mtrx[2][i]  =   satpos_rot_corr[2][i] - pos_ref_xyz[2];
        }
        double[] rho_predict =new double[nmbOfSatellites];//redundancy actually
        for(int i =0;i<nmbOfSatellites;i++)
        {
            rho_predict[i] = Math.sqrt(sat2usr_mtrx[0][i]*sat2usr_mtrx[0][i]+sat2usr_mtrx[1][i]*sat2usr_mtrx[1][i]+sat2usr_mtrx[2][i]*sat2usr_mtrx[2][i]);
        };
        double[] traveltime = new double[nmbOfSatellites];
        for(int i =0;i<nmbOfSatellites;i++)
        {
            traveltime[i]=rho_predict[i] / 299792458;
        }
        satpos_rot_corr = new double[3][nmbOfSatellites];
        for(int i =0;i<nmbOfSatellites;i++){                //Calculate the difference in signal propagation due to earth's rotation
            double omegatau   = 7.292115147e-5 * traveltime[i];
            double[][] R3 = new double[3][3];
            double []X_sat_rot = new double[3];
            double []X_sat = new double[3];
            R3[0][0]=Math.cos(omegatau);            R3[0][1]=Math.sin(omegatau);            R3[0][2]=0;
            R3[1][0]=-Math.sin(omegatau);            R3[1][1]=Math.cos(omegatau);            R3[1][2]=0;
            R3[2][0]=0;                                             R3[2][1]=0;                                           R3[2][2]=1;
            for(int j =0;j<3;j++) {
                X_sat[j]=satpos_ref_actv[j][i];
            }
            X_sat_rot = MatrixMultiply(R3,X_sat) ;
            for(int j =0;j<3;j++) {
                satpos_rot_corr[j][i] = X_sat_rot[j];
            }
        }
        for(int i =0;i<nmbOfSatellites;i++) {
            sat2usr_mtrx[0][i]  =   satpos_rot_corr[0][i] - pos_ref_xyz[0];
            sat2usr_mtrx[1][i]  =   satpos_rot_corr[1][i] - pos_ref_xyz[1];
            sat2usr_mtrx[2][i]  =   satpos_rot_corr[2][i] - pos_ref_xyz[2];
        }
        for(int i =0;i<nmbOfSatellites;i++)
        {
            rho_predict[i] = Math.sqrt(sat2usr_mtrx[0][i]*sat2usr_mtrx[0][i]+sat2usr_mtrx[1][i]*sat2usr_mtrx[1][i]+sat2usr_mtrx[2][i]*sat2usr_mtrx[2][i]);
        };
        obs_corr_actv = new double[nmbOfSatellites];

        obs_corr_actv = MinusArray(obs_actv,obs_ref_actv);

        double[]obs_corr_actv_DD = new double[nmbOfSatellites-1];
        for(int i =0 ;i< (nmbOfSatellites);i++){
            if (i < max_el){
                obs_corr_actv_DD[i] = obs_corr_actv[i] - obs_corr_actv[max_el];
            }else if(i > max_el){
                obs_corr_actv_DD[i-1] = obs_corr_actv[i] - obs_corr_actv[max_el];
            }
        }
        H = new double[nmbOfSatellites][3];
        H_DD = new double[nmbOfSatellites-1][3];
        for(int i=0;i<nmbOfSatellites;i++)
        {
            H[i][0] = -sat2usr_mtrx[0][i]/rho_predict[i];  //Ir = H
            H[i][1] = -sat2usr_mtrx[1][i]/rho_predict[i];
            H[i][2] = -sat2usr_mtrx[2][i]/rho_predict[i];
        }
        for(int i =1 ;i< (nmbOfSatellites);i++){
            if (i < max_el){
                for(int j =0;j<3;j++)
                {
                    H_DD[i][0] = H[i][0] - H[max_el][0];
                    H_DD[i][1] = H[i][1] - H[max_el][1];
                    H_DD[i][2] = H[i][2] - H[max_el][2];
                }
            }else if(i > max_el){
                for(int j =0;j<3;j++)
                {
                    H_DD[i-1][0] = H[i][0] - H[max_el][0];
                    H_DD[i-1][1] = H[i][1] - H[max_el][1];
                    H_DD[i-1][2] = H[i][2] - H[max_el][2];
                }
            }
        }

        deltaP_DD=new double[nmbOfSatellites];
        sat_vel_DD=new double[nmbOfSatellites];

        deltaP_DD = MinusArray(deltaP,deltaP_ref);
        double[][]temp = new double[3][nmbOfSatellites];
        temp =MinusArray(satvel_actv,satvel_ref_actv);
        for(int i =0;i<nmbOfSatellites;i++)
        {
            for(int j =0;j<3;j++)
            {
                sat_vel_DD[i]+=H[i][j]*   temp[j][i]   ;
            }
        }
        double[]obs_vel_DD = new double[nmbOfSatellites-1];
        for(int i =0 ;i< (nmbOfSatellites);i++)
        {
            if (i < max_el){
                obs_vel_DD[i] =deltaP_DD[i]-deltaP_DD[max_el]+sat_vel_DD[i] -sat_vel_DD[max_el];
            }else if(i > max_el){
                obs_vel_DD[i-1] = deltaP_DD[i]-deltaP_DD[max_el] +sat_vel_DD[i] -sat_vel_DD[max_el];
            }
        }
        double[] rhodot_predict ;
        double[] Rur ;
        double[] omc_rho =  new double[nmbOfSatellites-1];
        double[] omc_rhodot =  new double[nmbOfSatellites-1];
        rhodot_predict  = MatrixMultiply(H_DD,nxtState_vel);

        Rur = MatrixMultiply(H_DD,bur);
        omc_rho = MinusArray(obs_corr_actv_DD , Rur);//observation - prediction
        omc_rhodot =MinusArray( obs_vel_DD,rhodot_predict);

        double[]drho = new double[2*(nmbOfSatellites -1)];
        double[] nxtstate_ = new double[6];
        nxtstate_[0]  =KalFilt.stt_x[0];        nxtstate_[1] =KalFilt.stt_x[1];
        nxtstate_[2]  =KalFilt.stt_y[0];        nxtstate_[3]  =KalFilt.stt_y[1];
        nxtstate_[4]  =KalFilt.stt_z[0];        nxtstate_[5]  =KalFilt.stt_z[1];
        double[][] HK = new double[2*(nmbOfSatellites -1)][6];
        double[][] Rdiag = new double[2*(nmbOfSatellites -1)][2*(nmbOfSatellites -1)];

        for(int n = 0 ; n<nmbOfSatellites-1;n++)
        {
            drho[n*2] =omc_rho[n];
            drho[n*2+1] =omc_rhodot[n];
            HK[n*2][0]  =  H_DD[n][0];
            HK[n*2][2]  =  H_DD[n][1];
            HK[n*2][4]  =  H_DD[n][2];
            HK[n*2+1][1]  =  H_DD[n][0];
            HK[n*2+1][3]  =  H_DD[n][1];
            HK[n*2+1][5]  =  H_DD[n][2];
            double[] temp1 = new double [2];
            temp1 = EKF_R_Compute_new1( LLH[0], el[n], cn0_actv[n], Rv);

            Rdiag[n*2][n*2] = temp1[0];
            Rdiag[n*2+1][n*2+1] = temp1[1];
        }
        Pk = KalFilt.P;
        double[][] Kgain = new double[6][2*(nmbOfSatellites -1)];
        double[] dx = new double[6];
        Kgain =MatrixMultiply( MatrixMultiply(Pk,MatrixTranspose( HK)),
                getReverseMatrix(AddArray(MatrixMultiply(MatrixMultiply(HK,Pk),MatrixTranspose( HK)),Rdiag)));
        dx = MatrixMultiply(Kgain,drho);

        KalFilt.P = MatrixMultiply((  MinusArray(MatrixIdentify(6), MatrixMultiply(Kgain, HK)  )) ,Pk) ;
        double[] newState = new double[6];
        newState = AddArray(nxtstate_,dx);
        for(int i = 0;i<2;i++)
        {
            KalFilt.stt_x[i] = newState[i];
            KalFilt.stt_y[i]  = newState[i+2];
            KalFilt.stt_z[i]  = newState[i+4];
        }

        bur[0] = newState[0];
        bur[1] = newState[2];
        bur[2] = newState[4];
        //TODO: to define which parameters should be return

        // for(int i = 0;i<bur.length;i++){
        //                        System.out.println("i="+i+"  origin Function= "+ bur[i]);             
        // }
    }

    public double[] getBur(PvtCalculator pvc){
        double[] temp2 = new double[3];

            temp2[0] = pvc.stt_x[0];
            temp2[1] = pvc.stt_y[0];
            temp2[2] = pvc.stt_z[0];

            return temp2;
    }
    public double[] getBur(){
        // System.out.println(bur[0]);
            return bur;
    }


    public void Predict(PvtCalculator pvtCalculator)//TODO:timestamp is default
    {
        pvtCalculator.stt_x[0] = pvtCalculator.stt_x[0]+pvtCalculator.stt_x[1];
        pvtCalculator.stt_y[0] = pvtCalculator.stt_y[0]+pvtCalculator.stt_y[1];
        pvtCalculator.stt_z[0] = pvtCalculator.stt_z[0]+pvtCalculator.stt_z[1];
        double[][] A=MatrixIdentify(6);
        for(int i = 0 ; i<6;i++)
        {
            if(i%2==0)
            {
                A[i][i+1] = 1;
            }
        }
         pvtCalculator.P = AddArray(MatrixMultiply( MatrixMultiply(A,pvtCalculator.P), MatrixTranspose(A)), Qw);
    }

    public double[] CopyArray(double[] array1)
    {
        double[] array2 = new double[array1.length];
        for(int i =0;i<array1.length;i++)
        {
            array2[i] = array1 [i];
        }
        return array2;
    }

    public double[] AddArray(double[] array1,double[] array2)
    {
        double[] array3 = new double[array1.length];
        for(int i =0;i<array1.length;i++)
        {
            array3[i] = array1 [i]+ array2[i];
        }
        return array3;
    }
    public double[][] AddArray(double[][] array1,double[][] array2)
    {
        double[][] array3 = new double[array1.length ][array1[1].length];
        for(int i =0;i<array1.length;i++)
        {
            for(int j=0;j<array1[1].length;j++)
            {
                array3[i][j] = array1 [i][j]+array2[i][j];
            }

        }
        return array3;
    }

    public double[] MinusArray(double[] array1,double[] array2)
    {
        double[] array3 = new double[array1.length];
        for(int i =0;i<array1.length;i++)
        {
            array3[i] = array1 [i]- array2[i];
        }
        return array3;
    }
    public double[][] MinusArray(double[][] array1,double[][] array2)
    {
        double[][] array3 = new double[array1.length ][array1[1].length];
        for(int i =0;i<array1.length;i++)
        {
            for(int j=0;j<array1[1].length;j++)
            {
                array3[i][j] = array1 [i][j]- array2[i][j];
            }

        }
        return array3;
    }

    public double[][]  MatrixMultiply(double[][] array1,double[][] array2)
    {

        int n1 = array1.length;
        int m1 = array1[1].length;
        int n2 = array2.length;
        int m2 = array2[1].length;
        double[][] array3 = new double[n1 ][m2];
        if(m1 != n2)
        {
       }else
        {
            for(int i = 0;i<n1;i++)
            {
                for (int j = 0; j < m2; j++)
                {
                    for(int n = 0;n<m1;n++)
                    {
                        array3[i][j] += array1[i][n] * array2[n][j];
                    }

                }
            }
        }
        return array3;
    }

    public double[]  MatrixMultiply(double[][] array1,double[] array2)
    {

        int n1 = array1.length;
        int m1 = array1[1].length;
        int n2 = array2.length;
        int m2 = 1;
        double[] array3 = new double[n1 ];
        if(m1 != n2)
        {
       }else
        {
            for(int i = 0;i<n1;i++)
            {
                    for(int n = 0;n<m1;n++)
                    {
                        array3[i] += array1[i][n] * array2[n];
                    }

            }
        }
        return array3;
    }

    public double[][]  MatrixBlock(double[][] array1,int n)
    {

        int n1 = array1.length;
        int m1 = array1[0].length;
        int m2 = 1;
        double[][] array3 = new double[n1*n ][m1*n];

        for(int j = 0;j<n;j++)
        {
            for(int i =0;i< n1;i++)
            {
                for(int m =0;m<m1;m++)
                {
                    array3[j*n1+i][j*m1+m] = array1[i][m];
                }

            }
        }
        return array3;
    }

    public double[][]  MatrixTranspose(double[][] array1)
    {
        double[][] A = new double[array1[0].length][array1.length];
        for (int i = 0; i < array1.length; i++){
            for (int j = 0; j < array1[0].length; j++){
                A[j][i] = array1[i][j];
            }
        }
        return A;
    }
    public double[][]  MatrixIdentify(int n)
    {
        double[][] A = new double[n][n];
        for (int i = 0; i <n; i++)
                A[i][i] = 1;
        return A;
    }
    public int find(double[] array1,int tag)
    {
        if(tag == 1) //tag=1 找最大的MAX的position
        {
            double max;
            int position;
            max = array1[0];
            position = 0;
            for(int i = 1;i<array1.length;i++)
            {
                if (array1[i]>max)
                {
                    max = array1[i];
                    position = i;
                }
            }
        return position;
        }else{
            return 0 ;
        }
    }
    /** algebraic complement
     * @param h row
     * @param v column
     */
    public  double[][] getDY(double[][] data, int h, int v) {
        int H = data.length;
        int V = data[0].length;

        double[][] newData = new double[H - 1][V - 1];

        for (int i = 0; i < newData.length; i++) {

            if (i < h - 1) {
                for (int j = 0; j < newData[i].length; j++) {
                    if (j < v - 1) {
                        newData[i][j] = data[i][j];
                    } else {
                        newData[i][j] = data[i][j + 1];
                    }
                }
            } else {
                for (int j = 0; j < newData[i].length; j++) {
                    if (j < v - 1) {
                        newData[i][j] = data[i + 1][j];
                    } else {
                        newData[i][j] = data[i + 1][j + 1];
                    }
                }

            }
        }
        return newData;
    }
    /**
     * Magnitude of determinant
     */
    public  double getHL(double[][] data) {

        if (data.length == 2) {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        }
        if (data.length == 1) {
            return data[0][0];
        }
        double total = 0;
        int num = data.length;
        double[] nums = new double[num];

        for (int i = 0; i < num; i++) {
            if (i % 2 == 0) {
                nums[i] = data[0][i] * getHL(getDY(data, 1, i + 1));
            } else {
                nums[i] = -data[0][i] * getHL(getDY(data, 1, i + 1));
            }
        }

        for (int i = 0; i < num; i++) {
            total += nums[i];
        }
        return total;
    }
    /**
     * Get the reverse matrix
     */
    public  double[][] getReverseMatrix(double[][] data) {

        double m = getHL(data);
        if(m == 0){
        }
        double[][] newData = new double[data.length][data.length];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data.length; j++) {
                double num;
                if ((i + j) % 2 == 0) {
                    num = getHL(getDY(data, i + 1, j + 1));
                } else {
                    num = -getHL(getDY(data, i + 1, j + 1));
                }
                newData[i][j] = num / m;
            }
        }
        newData = MatrixTranspose(newData);
        return newData;
    }


    public double[]   ECEF2geo(double x,double y, double z)  // the height calculation is a bit difference from matlab but no effect.
    {

   // Convert from ECEF cartesian coordinates to  latitude, longitude and height.  WGS-84
        double x2 = x *x;
        double y2 = y *y;
        double z2 = z *z;

        double a = 6378137.0000 ;   //earth radius in meters
        double b = 6356752.3142;    // earth semiminor in meters
        double e = Math.sqrt (1-(b/a)*(b/a));
        double b2 = b*b;
        double e2 = e *e ;
        double ep = e*(a/b);
        double r = Math.sqrt(x2+y2);
        double r2 = r*r;
        double E2 = a *a - b *b;
        double F = 54*b2*z2;
        double G = r2 + (1-e2)*z2 - e2*E2;
        double c = (e2*e2*F*r2)/(G*G*G);
        double s =Math.pow(  ( 1 + c + Math.sqrt( c*c + 2*c) ),(1/3));
        double P = F / (3 * Math.pow((s+1/s+1),2)  * G*G);
        double Q = Math.sqrt(1+2*e2*e2*P);
        double ro = -(P*e2*r)/(1+Q) + Math.sqrt((a*a/2)*(1+1/Q) - (P*(1-e2)*z2)/(Q*(1+Q)) - P*r2/2);
        double tmp =Math.pow((r - e2*ro),2);
        double U = Math.sqrt( tmp + z2 );
        double V = Math.sqrt( tmp + (1-e2)*z2 );
        double zo = (b2*z)/(a*V);
        double height = U*( 1 - b2/(a*V) );
        double  lat = Math.atan( (z + ep*ep*zo)/r );
        double  temp = Math.atan(y/x);
        double [] LLH = new double[3];
        if(x>=0)
        {
              LLH[1]= temp;
        }else if((x < 0) & (y >= 0)){
              LLH[1] = PI+ temp;
        }else{
            LLH[1] =  temp - PI;
        }
        LLH[0] = lat/(PI/180);
        LLH[1]  = LLH[1] /(PI/180);
        LLH[2] = height;

        return LLH;
    }

    public double[] EKF_R_Compute_new1(double lat, double el_actv, double cn0_actv, double[] Rv0)
    {
        double Re = 6378.14; //earch radius, km
        double hI = 350;     //height of maximum electron density, km
        double d2r = PI/180;
        double Bn = 1;
        double BL = 10;
        double T  = 0.02;

        long L1= 1575420000;
        double B1= 1561.098e6;
        double cTc = C/1.023e6;
        double cf0 = C/L1;
        double sigma_v,sigma_t;
        if( lat<20){
            sigma_v = 9;
        }else if(lat <55){
            sigma_v = 4.5;
        }else{
            sigma_v = 6;
        }
        sigma_t = 0.12;
        double[] Rvec = new double[2];
        double sigma2_iono = sigma_v / (1 + Math.pow( Re*Math.sin((d2r*el_actv) )/ (Re+hI) ,2));
        double sigma2_trop = Math.pow(sigma_t,2) * Math.pow(1.001 ,2) / (0.002001 +     Math.pow( (Math.sin(d2r*el_actv)),2));

        double cn0_lin = Math.pow(10,(cn0_actv/10));
        double sigma2_psr = Math.pow(cTc,2)* (Bn/cn0_lin/2) * (1 + 2/T/cn0_lin) + Rv0[0];
        double sigma2_psrdot=Math.pow( (cf0/T/2/PI),2) * (4*BL/cn0_lin * (1 + 1/T/cn0_lin)) + Rv0[1];
//        if cn0_actv(2,n)  //this line origin purpose is to judge BDS or GPS
//        double cn0_los2mp = cn0_actv - cn0_actv;
//        double sigma2_mpratio = 10^(-1.1*(cn0_los2mp - 10)/10);
//    else
        double sigma2_mpratio = 0;
        Rvec[0] = sigma2_iono + sigma2_trop + sigma2_psr*(1+sigma2_mpratio);
        Rvec[1]   = sigma2_psrdot;
        return Rvec;
    }


    public static void main(String[] args) {
        DDKalmanFilterManager A = new DDKalmanFilterManager();
        int nmbOfSatellites =5;
        int[] activeChannel   = new int[nmbOfSatellites]; // activeChannel array saves the prn of every GPS Satellite 

        double[][] satpos  =new double[6][32] ;// 6 rows reserve the sat position_xyz and velocity_xyz; 32 volumns reserve all the 32 prn;
        double[][] satpos_ref = new double[6][32] ;

        double[] obs = new double[32];//every 32 GPS satellite pseudorange, invalid satellite ele is zer
        double[] obs_ref = new double[32];

        double[] doppSmooth =new double[32]; // //every 32 GPS satellite doppler value, invalid satellite ele is zer
        double[] doppSmooth_ref =new double[32];

        double[] elfore = new double[32]; //every 32 GPS satellite elevation.invalid satellite ele is zero
        double[] azfore = new double[32]; //every 32 GPS satellite elevation.invalid satellite ele is zero
        double[]cn0 = new double[32]; //every 32 GPS satellite elevation.invalid satellite ele is zero

        PvtCalculator pvc = new PvtCalculator();
        // 以下取matlab第九个loop的值 在kalmanPos_LOG里面
        pvc.stt_x[0] =142.666079794193;
        pvc.stt_x[1] =-0.918552094629031;
        pvc.stt_y[0] =147.342708506158;
        pvc.stt_y[1] =-0.343336683501094;
        pvc.stt_z[0] =-113.957337514906;
        pvc.stt_z[1] =0.408781147882364;
        pvc.P[0][0] = 3.57139274266416;pvc.P[0][1] = 0.357833891077047;
        pvc.P[1][0] = 0.357833891077047;pvc.P[1][1] = 0.583693644032989;

        pvc.P[2][2] = 4.03560566100503;pvc.P[2][3] = 0.365795584604125;
        pvc.P[3][3] = 0.365795584604125;pvc.P[3][4] = 0.584678501596958;

        pvc.P[4][4] = 6.87951773763442;pvc.P[4][5] = 0.522202696919255;
        pvc.P[5][4] = 0.522202696919255;pvc.P[5][5] = 0.677924484233650;

        activeChannel[0] = 10;
        activeChannel[1] = 12;
        activeChannel[2] = 14;
        activeChannel[3] = 20;
        activeChannel[4] = 24;

        elfore[activeChannel[0]-1] = 73.7625208523695;
        elfore[activeChannel[1]-1] = 39.2206528131214;
        elfore[activeChannel[2]-1] = 31.1380302421402;
        elfore[activeChannel[3]-1] = 50.3160322717088;
        elfore[activeChannel[4]-1] = 17.4308330569214;

        azfore[activeChannel[0]-1] = 232.359978068975;
        azfore[activeChannel[1]-1] = 58.7483580945080;
        azfore[activeChannel[2]-1] = 307.925855479008;
        azfore[activeChannel[3]-1] = 179.578154780576;
        azfore[activeChannel[4]-1] = 52.3870389691604;

        doppSmooth[activeChannel[0]-1] = 43.4691642646761;
        doppSmooth[activeChannel[1]-1] = 247.210890919701;
        doppSmooth[activeChannel[2]-1] = -447.471385124554;
        doppSmooth[activeChannel[3]-1] = 441.172093673909;
        doppSmooth[activeChannel[4]-1] = 521.384302044530;

        doppSmooth_ref[activeChannel[0]-1] = 89.8802687108149;
        doppSmooth_ref[activeChannel[1]-1] = 294.941111105679;
        doppSmooth_ref[activeChannel[2]-1] = -401.642389023505;
        doppSmooth_ref[activeChannel[3]-1] = 487.377110072407;
        doppSmooth_ref[activeChannel[4]-1] = 569.584737895992;

        obs[activeChannel[0]-1] = 20253750.6200000;
        obs[activeChannel[1]-1] = 21805007.2080000;
        obs[activeChannel[2]-1] = 23018180.0460000;
        obs[activeChannel[3]-1] = 21106270.1330000;
        obs[activeChannel[4]-1] = 23978486.3390000;

        obs_ref[activeChannel[0]-1] = 20254118.8910000;
        obs_ref[activeChannel[1]-1] = 21805141.2810000;
        obs_ref[activeChannel[2]-1] = 23018560.1560000;
        obs_ref[activeChannel[3]-1] = 21106659.7730000;
        obs_ref[activeChannel[4]-1] = 23978574.8830000;


        satpos[0][activeChannel[0]-1] =-8661379.76080410;
        satpos[1][activeChannel[0]-1] =22760331.6193311;
        satpos[2][activeChannel[0]-1] =10436610.4446874;
        satpos[3][activeChannel[0]-1] =-1021.21588229509;
        satpos[4][activeChannel[0]-1] =882.000655558191;
        satpos[5][activeChannel[0]-1] =-2798.56859505564;

        satpos[0][activeChannel[1]-1] =-19009036.2901217;
        satpos[1][activeChannel[1]-1] =3287133.31467227;
        satpos[2][activeChannel[1]-1] =17986157.6381967;
        satpos[3][activeChannel[1]-1] =1458.68218530203;
        satpos[4][activeChannel[1]-1] =-1813.04372578464;
        satpos[5][activeChannel[1]-1] =1870.14878984547;

        satpos[0][activeChannel[2]-1] =8349700.40466628;
        satpos[1][activeChannel[2]-1] =16074829.7704544;
        satpos[2][activeChannel[2]-1] =19787086.8452291;
        satpos[3][activeChannel[2]-1] =-2404.46016697464;
        satpos[4][activeChannel[2]-1] =-436.853256066700;
        satpos[5][activeChannel[2]-1] =1373.63144723353;

        satpos[0][activeChannel[3]-1] =-13924621.2610703;
        satpos[1][activeChannel[3]-1] =22551058.3861175;
        satpos[2][activeChannel[3]-1] =163336.039620139;
        satpos[3][activeChannel[3]-1] =-349.822180974318;
        satpos[4][activeChannel[3]-1] =-176.941770291960;
        satpos[5][activeChannel[3]-1] =-3110.16856309030;

        satpos[0][activeChannel[4]-1] =-17711421.2111442;
        satpos[1][activeChannel[4]-1] =-5730945.86273570;
        satpos[2][activeChannel[4]-1] =18930294.2474725;
        satpos[3][activeChannel[4]-1] =-857.567696953741;
        satpos[4][activeChannel[4]-1] =-2274.54740874787;
        satpos[5][activeChannel[4]-1] =-1451.98201510106;

        satpos_ref[0][activeChannel[0]-1] =-8661380.05028981;
        satpos_ref[1][activeChannel[0]-1] =22760331.8693532;
        satpos_ref[2][activeChannel[0]-1] =10436609.6513727;
        satpos_ref[3][activeChannel[0]-1] =-1021.21580626605;
        satpos_ref[4][activeChannel[0]-1] =882.000593807528;
        satpos_ref[5][activeChannel[0]-1] =-2798.56865847453;

        satpos_ref[0][activeChannel[1]-1] =-19009035.8754872;
        satpos_ref[1][activeChannel[1]-1] =3287132.79930954;
        satpos_ref[2][activeChannel[1]-1] =17986158.1697917;
        satpos_ref[3][activeChannel[1]-1] =1458.68219877868;
        satpos_ref[4][activeChannel[1]-1] =-1813.04380158420;
        satpos_ref[5][activeChannel[1]-1] =1870.14867876847;

        satpos_ref[0][activeChannel[2]-1] =8349699.72316517;
        satpos_ref[1][activeChannel[2]-1] =16074829.6466362;
        satpos_ref[2][activeChannel[2]-1] =19787087.2345603;
        satpos_ref[3][activeChannel[2]-1] =-2404.46022130307;
        satpos_ref[4][activeChannel[2]-1] =-436.853226503515;
        satpos_ref[5][activeChannel[2]-1] =1373.63133143499;

        satpos_ref[0][activeChannel[3]-1] =-13924621.3602100;
        satpos_ref[1][activeChannel[3]-1] =22551058.3359722;
        satpos_ref[2][activeChannel[3]-1] =163335.158197448;
        satpos_ref[3][activeChannel[3]-1] =-349.822124779304;
        satpos_ref[4][activeChannel[3]-1] =-176.941858685987;
        satpos_ref[5][activeChannel[3]-1] =-3110.16856408179;

        satpos_ref[0][activeChannel[4]-1] =-17711421.4550405;
        satpos_ref[1][activeChannel[4]-1] =-5730946.50962752;
        satpos_ref[2][activeChannel[4]-1] =18930293.8345220;
        satpos_ref[3][activeChannel[4]-1] =-857.567710813928;
        satpos_ref[4][activeChannel[4]-1] =-2274.54734713549;
        satpos_ref[5][activeChannel[4]-1] =-1451.98212977508;

        cn0[activeChannel[0]-1] = 38.6750000000000;
        cn0[activeChannel[1]-1] = 35.6590000000000;
        cn0[activeChannel[2]-1] = 35.5060000000000;
        cn0[activeChannel[3]-1] = 34.0960000000000;
        cn0[activeChannel[4]-1] = 22.1230000000000;
        double[] e  = new double[3];
        e = A.getBur(pvc);
        A.Predict(pvc);
        A.Update(activeChannel,satpos, satpos_ref, obs,obs_ref, pvc,  doppSmooth,doppSmooth_ref, elfore,azfore,cn0);
        e = A.getBur();
        for(int i=0; i<3;i++){
            System.out.println("bur="+e[i]);
        }
    }



 
}
