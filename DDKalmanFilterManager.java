package com.dzd.gnss2.Method;

/**
 * Created by HP on 2019/5/13.
 */
import android.location.GnssClock;
import android.location.GnssMeasurement;
import android.location.GnssStatus;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.DiagonalMatrix;


import com.dzd.gnss2.Data.NavData.GpsNavData;
import com.dzd.gnss2.Data.NavData.GpsNavL1CA;
import com.dzd.gnss2.Data.NavData.NavData;
import com.dzd.gnss2.Data.SatelliteData.SatelliteData;
import com.dzd.gnss2.Data.SatelliteData.SatelliteMeasurementStatus;
import com.dzd.gnss2.Other.Position;
import com.dzd.gnss2.Other.PositionConfig;

import java.util.ArrayList;
import java.util.Arrays;

import Jama.Matrix;

import static com.dzd.gnss2.Other.Constant.F;
import static com.dzd.gnss2.Other.Constant.OMEGA_dot_e;
import static com.dzd.gnss2.Other.Constant.mu;
import static java.lang.Math.floor;

public class DDKalmanFilterManager {

    private double nxtState_pos[] = new double[3];
    private double nxtState_vel[] = new double[3];
    private double Ac[][] = {{1,0},{0,1}};

    private static int Sp = 1;
    private static double Sv=0.5;
    private static double Pxyz0[] = {30,3,30,3,30,3} ;
    private double Qp0[][] = {{Sp+Sv/3, Sv/2},{Sv/2,Sv}};
    private double Pk[][] = new double[6][6];
    private double Qw[][] = new double[6][6];
    RealMatrix matrix_nxtState_pos = new Array2DRowRealMatrix(nxtState_pos);
    RealMatrix matrix_nxtState_vel = new Array2DRowRealMatrix(nxtState_vel);
    RealMatrix matrix_Qp = new Array2DRowRealMatrix(Qp0);
    RealMatrix matrix_Pk = new Array2DRowRealMatrix(Pk);
    RealMatrix matrix_Qw = new Array2DRowRealMatrix(Qw);

    private double  satpos_ref_actv[][];
    private double  satvel_ref_actv[][] ;
    private double  satpos_rot_corr[][]; //storing the sat positions after earth rotation corrections
    private double obs_corr_actv[] ;

    private double transmitTime_actv[] ;
    private double satClkCorr_actv[][] ;
    private double satClkCorr_ref_actv[][] ;
    private double cn0_actv[] ;
    private double  az[] ;
    private double  el[] ;
    private double trop[];
    private double iono[] ;
    private double DOP[] ;
    private double H[][] ;

    public  static double[] pos_ref_xyz= {-2853445.926, 4667466.476, 3268291.272};
    public  class  PvtCalculator{
        public double[] stt_x  = new double[2];
        public double[] stt_y = new double[2];
        public double[] stt_z  = new double[2];

    };

    DDKalmanFilterManager()
    {
        matrix_Pk =new DiagonalMatrix(Pxyz0);
        matrix_Qw =new DiagonalMatrix(Qp0);//初始化有误
    }

    public void Update(int[] activeChannel,double[][] satpos,double[] satpos_ref,double[] obs,double[] obs_ref, double[] transmitTime, PvtCalculator pvtCalculator, PvtCalculator pvtCalculator_ref,boolean pvtForecast_Succ, double[] elfore, double[] azfore,double[]cn0)
    {
        Initialize(activeChannel,satpos,] satpos_ref,  obs,  obs_ref,   transmitTime,    pvtCalculator,   pvtCalculator_ref, pvtForecast_Succ,   elfore,   azfore, cn0)
        int  nmbOfSatellites = activeChannel.length;
        double max = Arrays.stream(el).max().getAsDouble();
        PvtCalculator KalFilt =  pvtCalculator;
        nxtState_pos = new double[]{KalFilt.stt_x[0], KalFilt.stt_y[0], KalFilt.stt_z[0]};
        nxtState_vel = new double[]{KalFilt.stt_x[1], KalFilt.stt_y[1], KalFilt.stt_z[1]};

        double[] bur = nxtState_pos ;
        for (int i = 0;i<3;i++){
            nxtState_pos[i] = nxtState_pos[i] + pos_ref_xyz[i];
        }
        [ Lat, Lon, Hight ] = cart2geo( nxtState_pos[0], nxtState_pos[1], nxtState_pos[2], 5 );

        double [][] sat2usr_mtrx = new double[3][nmbOfSatellites]
        for(int i =0;i<nmbOfSatellites;i++) {
            sat2usr_mtrx[0][i]  =   satpos_rot_corr[0][i] - pos_ref_xyz[0];
            sat2usr_mtrx[1][i]  =   satpos_rot_corr[1][i] - pos_ref_xyz[1];
            sat2usr_mtrx[2][i]  =   satpos_rot_corr[2][i] - pos_ref_xyz[2];
        }
        double[] rho_predict =new double[nmbOfSatellites];
        for(int i =0;i<nmbOfSatellites;i++)
        {
            rho_predict[i] = Math.sqrt(sat2usr_mtrx[0][i]*sat2usr_mtrx[0][i]+sat2usr_mtrx[1][i]*sat2usr_mtrx[1][i]+sat2usr_mtrx[2][i]*sat2usr_mtrx[2][i]);
        };
        double[] traveltime = new double[nmbOfSatellites];
        for(int i =0;i<nmbOfSatellites;i++){
            traveltime[i]=rho_predict[i] / 299792458;
        }
        satpos_rot_corr = new double[3][nmbOfSatellites];
        for(int i =0;i<nmbOfSatellites;i++){
            for(int j =0;j<3;j++) {
                satpos_rot_corr[j][i] = e_r_corr(traveltime[i], satpos_ref_actv(j,i)); //解算地球自转的差异
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
        for(int i =0;i<nmbOfSatellites;i++)
        {
            obs_corr_actv[i] =obs[i] - obs_ref[i];
       };

        for(int i =1 ;i< (nmbOfSatellites);i++){
            if i ~ = max_el
            obs_corr_actv_DD(i) = obs_corr_actv(i) - obs_corr_actv(max_el);
            end
                    end
            obs_corr_actv_DD(max_el) = [];
        }

        // TODO:伪距差分


        // TODO:多普勒差分
     /*   gpsDoppSmooth = pvtCalculator.GPS.doppSmooth;
        gpsDoppSmooth_ref = pvtCalculator_ref.GPS.doppSmooth;
        deltaP = zeros(1, nmbOfSatellites);
        deltaP_ref = zeros(1, nmbOfSatellites);
        obs_vel = zeros(1, nmbOfSatellites);
        C = 299792458;
% L1= 1575420000;
% wavelengthL1 = C/L1;

        for n=1:nmbOfSatellites
        if gpsDoppSmooth(activeChannel(2,n),4) > 5e10
        deltaP(n) = (gpsDoppSmooth(activeChannel(2,n),1) - gpsDoppSmooth(activeChannel(2,n),2)) / pvtCalculator.pvtT; %积分多普勒一秒的变化量（m）由于前面用的dopper，所以此行无效
    else
%         deltaP(n) = -wavelengthL1*gpsDoppSmooth(activeChannel(2,n),3);
        deltaP(n) = gpsDoppSmooth(activeChannel(2,n),3);
        end

        if gpsDoppSmooth_ref (activeChannel(2,n),4) > 5e10
        deltaP_ref(n) = (gpsDoppSmooth_ref (activeChannel(2,n),1) - gpsDoppSmooth_ref(activeChannel(2,n),2)) / pvtCalculator.pvtT; %积分多普勒一秒的变化量（m）
    else
%         deltaP(n) = -wavelengthL1*gpsDoppSmooth(activeChannel(2,n),3);
        deltaP_ref(n) = gpsDoppSmooth_ref(activeChannel(2,n),3);
        end
        deltaP_DD(n) = deltaP(n)-deltaP_ref(n);
        sat_vel_DD(n)=H(n, 1:3)*(satvel_actv(:,n)-satvel_ref_actv(:,n));
%     obs_vel(n) = deltaP(n) - deltaP_ref(n) +H(n, 1:3)*(satvel_actv(:,n)-satvel_ref_actv(:,n));%本次根据定位方程得到的多普勒值（在单位矢量方向上的多普值单位m/s)
        end
        for i= 1 : nmbOfSatellites
        if i ~=max_el
        obs_vel_DD(i)  = deltaP_DD(i)-deltaP_DD(max_el) +sat_vel_DD(i) -sat_vel_DD(max_el);
        end
                end
        obs_vel_DD(max_el)= [];
        */
        // TODO:计算H矩阵
        /*
        drho = zeros(5*(nmbOfSatellites-1), 1); % [2*nmbOfSatellites x 1]
nxtstate_ = [KalFilt.stt_x; KalFilt.stt_y; KalFilt.stt_z; KalFilt.stt_dtf(:, 2)]; % vector [8x1]
HK = zeros(5*(nmbOfSatellites -1), 12);
Rdiag = zeros(5*(nmbOfSatellites -1), 1);
for n=1:(nmbOfSatellites-1)
    drho(n*5-4) = omc_rho(n);
    drho(n*5-3) = omc_rhodot(n);
    drho(n*5-2)   =  omc_acc(1);
        drho(n*5-1)   =  omc_acc(2);
            drho(n*5)   =  omc_acc(3);
    HK(n*5-4:n*5-2, :) = [H_DD(n,1)*eye(3), H_DD(n,2)*eye(3), H_DD(n,3)*eye(3), H_DD(n,4)*eye(3)];    %
    HK(n*5-2, :) = [0,0,1,0,0,0,0,0,0,0,0,0];
        HK(n*5-1, :) = [0,0,0,0,0,1,0,0,0,0,0,0];
            HK(n*5, :) = [0,0,0,0,0,0,0,0,1,0,0,0];
    Rdiag(n*5-4:n*5-3) = EKF_R_Compute_new1('GPS_L1CA', Lat, el(n), cn0_actv(:,n), KalFilt.Rv);
    Rdiag(n*5-2) = 0.1;
    Rdiag(n*5-1) = 0.1;
    Rdiag(n*5) = 0.1;
end
Pk_ = KalFilt.P; % Pk_ already computed in the function pvt_forecast_filt().
R = diag(Rdiag);

%------- 8th, performing the Kalman updating -------
Kgain = Pk_ * HK.' / (HK * Pk_ * HK.' + R);
dx = Kgain * drho;
KalFilt.P = (eye(12) - Kgain * HK) * Pk_;
KalFilt.P = (KalFilt.P + (KalFilt.P).')/2;    %滤波数值计算
newState = nxtstate_ + dx;

KalFilt.stt_x = newState(1:3);
KalFilt.stt_y = newState(4:6);
KalFilt.stt_z = newState(7:9);
KalFilt.stt_dtf(:, 2) = newState(10:12);

pos_xyz = [newState(1),newState(4),newState(7)]';
         */
    }

    public void Predict()
    {

    }

    public void Initialize(int[] activeChannel,double[][] satpos,double[][] satpos_ref,double[] obs,double[] obs_ref, double[] transmitTime, double[]  pvtCalculator, double[] pvtCalculator_ref,boolean pvtForecast_Succ, double[] elfore, double[] azfore,double[]cn0)
    {
        int  nmbOfSatellites = activeChannel.length;
        satpos_ref_actv =  new double[3][nmbOfSatellites];
        satvel_ref_actv=  new double[3][nmbOfSatellites];
        satpos_rot_corr =  new double[3][nmbOfSatellites]; //storing the sat positions after earth rotation corrections
        obs_actv = new double[nmbOfSatellites];
        obs_ref_actv = new double[nmbOfSatellites];
        cn0_actv= new double[nmbOfSatellites];
        az = new double[nmbOfSatellites];
        el = new double[nmbOfSatellites];;
        //       double DOP= new double[5];
        H= new double[4][nmbOfSatellites];

        for(int i =0;i<=nmbOfSatellites;i++)
        {
           // System.arraycopy(satpos_ref[i], 0, aa[i], 0, a[0].length);
            for(int j = 0;j<3;j++)
            {
                satpos_ref_actv[j][i] = satpos_ref[j][activeChannel[i]];
                satvel_ref_actv[j][i]= satpos_ref[j+3][activeChannel[i]];
            }
            obs_actv[i] = obs[activeChannel[i]];
            obs_ref_actv[i]=obs_ref[activeChannel[i]];
            cn0_actv[i] = cn0[ activeChannel[i]];
            el[i]   = elfore[activeChannel[i]];
            az[i]   = azfore[activeChannel[i]];
        }



    }

    public double cart2geo(double x, double y,double z){

    }






}
