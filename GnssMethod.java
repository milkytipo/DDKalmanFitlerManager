package com.dzd.gnss2.Method;

import android.location.GnssClock;
import android.location.GnssMeasurement;
import android.location.GnssStatus;


import com.dzd.gnss2.Data.NavData.GpsNavData;
import com.dzd.gnss2.Data.NavData.GpsNavL1CA;
import com.dzd.gnss2.Data.NavData.NavData;
import com.dzd.gnss2.Data.SatelliteData.SatelliteData;
import com.dzd.gnss2.Data.SatelliteData.SatelliteMeasurementStatus;
import com.dzd.gnss2.Other.Position;
import com.dzd.gnss2.Other.PositionConfig;

import java.util.ArrayList;
import Jama.Matrix;

import static com.dzd.gnss2.Other.Constant.F;
import static com.dzd.gnss2.Other.Constant.OMEGA_dot_e;
import static com.dzd.gnss2.Other.Constant.mu;
import static java.lang.Math.floor;

public class GnssMethod {

    /**
     *计算卫星位置 如果GPS存在L5CNAV电文则用L5CNAV计算位置坐标
     */
    static public Position calculateSvPos(GnssClock gnssClock, GnssMeasurement gnssMeasurement, NavData navData, SatelliteMeasurementStatus satelliteMeasurementStatus){
        //如果不能计算卫星位置 返回一个非法Position对象
        if(!satelliteMeasurementStatus.isCalculateSvPos()) return new Position();

        //根据不同的导航电文来计算 (因为可能算法会不一样，保留这种可能性，代码会冗长一点)
        if(satelliteMeasurementStatus.isGpsL5CNAVValid()){ //优先L5 CNAV

        }else if(satelliteMeasurementStatus.isGpsL1CAValid()){
            GpsNavL1CA gpsNavL1CA = ((GpsNavData) navData).getGpsNavL1CA();

            //计算E_k以及t_k
            double A = Math.pow(gpsNavL1CA.rootA,2);
            double n_0 = Math.pow(mu / Math.pow(A,3),0.5);
            double t = gnssMeasurement.getReceivedSvTimeNanos() / 1e9;
            double t_k = t - gpsNavL1CA.t_oe;
            if(t_k > 302400){
                t_k -= 604800;
            }else if(t_k < -302400){
                t_k += 604800;
            }
            double n = n_0 + gpsNavL1CA.delta_n;
            double M_k = gpsNavL1CA.M_0 + n * t_k;
            double e = gpsNavL1CA.e;
            double E_k = keplerSolve(e,M_k,0.00000000001);

            //计算钟差
            double t_sv = gnssMeasurement.getReceivedSvTimeNanos() / 1e9;
            double dif_tsv_toc = t_sv - gpsNavL1CA.t_oc;
            if(dif_tsv_toc > 302400){
                dif_tsv_toc -= 604800;
            }else if(dif_tsv_toc < -302400){
                dif_tsv_toc += 604800;
            }
            double delta_t_sv = gpsNavL1CA.a_f0 + gpsNavL1CA.a_f1 * dif_tsv_toc + gpsNavL1CA.a_f2 * dif_tsv_toc * dif_tsv_toc;
            delta_t_sv += F * gpsNavL1CA.e * gpsNavL1CA.rootA * Math.sin(E_k);
            //T_GD放到下一行代码中，不加在delta_t_sv中，和Correction类中的描述保持一致

            //加入钟差再次计算E_k
            t = gnssMeasurement.getReceivedSvTimeNanos() / 1e9 - delta_t_sv + gpsNavL1CA.T_GD;
            t_k = t - gpsNavL1CA.t_oe;
            if(t_k > 302400){
                t_k -= 604800;
            }else if(t_k < -302400){
                t_k += 604800;
            }
            n = n_0 + gpsNavL1CA.delta_n;
            M_k = gpsNavL1CA.M_0 + n * t_k;
            e = gpsNavL1CA.e;
            E_k = keplerSolve(e,M_k,0.00000000001);

            //计算卫星位置
            double sin_v_k = (Math.sqrt(1-Math.pow(gpsNavL1CA.e,2)) * Math.sin(E_k)) / (1 - gpsNavL1CA.e * Math.cos(E_k));
            double cos_v_k = (Math.cos(E_k) - gpsNavL1CA.e) / (1 - gpsNavL1CA.e * Math.cos(E_k));
            double v_k = Math.atan2(sin_v_k,cos_v_k);
            double PHI_k = v_k + gpsNavL1CA.omega;
            double delta_u_k = gpsNavL1CA.C_us * Math.sin(2 * PHI_k) + gpsNavL1CA.C_uc * Math.cos(2 * PHI_k);
            double delta_r_k = gpsNavL1CA.C_rs * Math.sin(2 * PHI_k) + gpsNavL1CA.C_rc * Math.cos(2 * PHI_k);
            double delta_i_k = gpsNavL1CA.C_is * Math.sin(2 * PHI_k) + gpsNavL1CA.C_ic * Math.cos(2 * PHI_k);
            double u_k = PHI_k + delta_u_k;
            double r_k = Math.pow(gpsNavL1CA.rootA,2) * (1 - gpsNavL1CA.e * Math.cos(E_k)) + delta_r_k;
            double i_k = gpsNavL1CA.i_0 + delta_i_k + gpsNavL1CA.IDOT * t_k;
            double x_k_aps = r_k * Math.cos(u_k);
            double y_k_aps = r_k * Math.sin(u_k);
            double OMEGA_k = gpsNavL1CA.OMEGA_0 + gpsNavL1CA.OMEGA_dot * t_k - OMEGA_dot_e * ((gnssClock.getTimeNanos() - gnssClock.getBiasNanos() - gnssClock.getFullBiasNanos()) % 604800000000000d + gnssMeasurement.getTimeOffsetNanos()) / 1e9;
            Position svPosition = new Position(x_k_aps * Math.cos(OMEGA_k) - y_k_aps * Math.cos(i_k) * Math.sin(OMEGA_k),
                    x_k_aps * Math.sin(OMEGA_k) + y_k_aps * Math.cos(i_k) * Math.cos(OMEGA_k),
                    y_k_aps * Math.sin(i_k));
            return svPosition;
        }else if(satelliteMeasurementStatus.isBeidouB1Valid()){

        }else {

        }
        return new Position();
    }

    /**
     * 计算用户位置
     */
    static public Position calculateUserPos(ArrayList<SatelliteData> satelliteDataArrayList, PositionConfig positionConfig){
        //选择用于定位的卫星测量量
        ArrayList<SatelliteData> selectedSatelliteDataArraylist = new ArrayList<SatelliteData>();
        for (int index = 0;index < satelliteDataArrayList.size();++index){
            boolean selectFlag = true;
            SatelliteData satelliteData = satelliteDataArrayList.get(index);
            //选择相应的卫星系统
            switch (positionConfig.constellation){
                case 0://GPS
                    if(satelliteData.gnssMeasurement.getConstellationType() != GnssStatus.CONSTELLATION_GPS)
                        selectFlag = false;
                    break;
                case 1://BEIDOU
                    if(satelliteData.gnssMeasurement.getConstellationType() != GnssStatus.CONSTELLATION_BEIDOU)
                        selectFlag = false;
                    break;
                case 2://GPS+BEIDOU
                    if(satelliteData.gnssMeasurement.getConstellationType() != GnssStatus.CONSTELLATION_BEIDOU
                            && satelliteData.gnssMeasurement.getConstellationType() != GnssStatus.CONSTELLATION_GPS)
                        selectFlag = false;
                    break;
                default:
                    selectFlag = false;
                    break;
            }


            //选择相应频段
            if(!positionConfig.dualFrequency && satelliteData.gnssMeasurement.hasCarrierFrequencyHz()){
                //只用GPS L1及 BEIDOU B1
                if(satelliteData.gnssMeasurement.getConstellationType() == GnssStatus.CONSTELLATION_GPS &&
                        Math.abs(satelliteData.gnssMeasurement.getCarrierFrequencyHz() - 1575.42e6) > 1e6)
                    selectFlag = false;
                if(satelliteData.gnssMeasurement.getConstellationType() == GnssStatus.CONSTELLATION_BEIDOU &&
                        Math.abs(satelliteData.gnssMeasurement.getCarrierFrequencyHz() - 1561.098e6) > 1e6)
                    selectFlag = false;
            }

            //选择pesudorange和卫星坐标都计算好了的卫星
            if(!satelliteData.satelliteMeasurementStatus.isCalculatePseudorange()
                    || !satelliteData.satelliteMeasurementStatus.isCalculateSvPos()){
                selectFlag = false;
            }

            if(selectFlag) selectedSatelliteDataArraylist.add(satelliteData);
        }

        //少于4个卫星 无法定位
        if(selectedSatelliteDataArraylist.size()<4) return new Position();

        //正式定位
        Position roughPosition,userPosition;
        try {
            //获取粗略定位结果
            roughPosition = getRoughPosition(selectedSatelliteDataArraylist);
            //精确定位(泰勒迭代法,以初略定位结果作为初值迭代)
            userPosition = getUserPosition(selectedSatelliteDataArraylist,roughPosition);
        }catch (RuntimeException e){
            roughPosition = new Position();
            userPosition = new Position();
        }
        return userPosition;
    }


    /**
     * 求粗略的position
     */
    static public Position getRoughPosition(ArrayList<SatelliteData> satelliteDataArrayList) {
        int sv_num = satelliteDataArrayList.size();
        Matrix[] a_i = new Matrix[sv_num];
        Matrix A = new Matrix(sv_num,4);
        Matrix i_0 = new Matrix(sv_num,1,1);
        Matrix r = new Matrix(sv_num,1);
        for (int sv_index = 0;sv_index < sv_num;++sv_index) {
            a_i[sv_index] = new Matrix(4,1);
            a_i[sv_index].set(0,0,satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_x);
            a_i[sv_index].set(1,0,satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_y);
            a_i[sv_index].set(2,0,satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_z);
            a_i[sv_index].set(3,0,satelliteDataArrayList.get(sv_index).getCorrectedPseudorange());
            A.setMatrix(sv_index,sv_index,0,3,a_i[sv_index].transpose());
            r.set(sv_index,0,0.5 * calculateMinkowskiFunc(a_i[sv_index],a_i[sv_index]));
        }
        Matrix B = A.inverse();
        Matrix u = B.times(i_0);
        Matrix v = B.times(r);
        double E = calculateMinkowskiFunc(u,u);
        double F = calculateMinkowskiFunc(u,v) - 1;
        double G = calculateMinkowskiFunc(v,v);
        if(4 * Math.pow(F,2) - 4 * E * G >= 0){
            double lamda_1 = (-2*F + Math.sqrt(4 * Math.pow(F,2) - 4 * E * G)) / (2 * E);
            double lamda_2 = (-2*F - Math.sqrt(4 * Math.pow(F,2) - 4 * E * G)) / (2 * E);
            Matrix y1 = (u.times(lamda_1)).plus(v);
            Matrix y2 = (u.times(lamda_2)).plus(v);
            double delta1 = 0;
            double delta2 = 0;
            double y[];
            for(int sv_index = 0;sv_index < sv_num;++sv_index) {
                delta1 += Math.abs(Math.sqrt(Math.pow(y1.get(0, 0) - satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_x,2) +
                        Math.pow(y1.get(1, 0) - satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_y,2) +
                        Math.pow(y1.get(2, 0) - satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_z,2)) -
                        y1.get(3,0) - satelliteDataArrayList.get(sv_index).getCorrectedPseudorange());
                delta2 += Math.abs(Math.sqrt(Math.pow(y2.get(0, 0) - satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_x,2) +
                        Math.pow(y2.get(1, 0) - satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_y,2) +
                        Math.pow(y2.get(2, 0) - satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_z,2)) -
                        y2.get(3,0) - satelliteDataArrayList.get(sv_index).getCorrectedPseudorange());
            }
            if(delta1 < delta2){
                y = y1.getRowPackedCopy();
                y[3] = -y[3];
            }else{
                y = y2.getRowPackedCopy();
                y[3] = -y[3];
            }
            Position roughPosition = new Position(y[0],y[1],y[2]);
            roughPosition.setClockBias(y[3]);
            return roughPosition;

        }else{
            return new Position();
        }
    }
    static private double calculateMinkowskiFunc(Matrix a,Matrix b){
        double result = a.get(0,0) * b.get(0,0) +
                a.get(1,0) * b.get(1,0) +
                a.get(2,0) * b.get(2,0) -
                a.get(3,0) * b.get(3,0);
        return result;
    }

    /**
     * 泰勒迭代求精确的定位结果
     */
    static private Position getUserPosition(ArrayList<SatelliteData> satelliteDataArrayList,Position roughPosition){
        if(!roughPosition.validPosition) return new Position();

        //TODO 加入权重
//        double[][] R = new double[satelliteDataArrayList.size()][satelliteDataArrayList.size()];
//        for(int sv_index = 0;sv_index < satelliteDataArrayList.size();++sv_index){
//            if(satelliteDataArrayList..URA != Double.MAX_VALUE) {
//                R[sv_index][sv_index] = Math.pow(navData[sv_index].URA, 2) + Math.pow(obsData.GPS_ReceivedSvTimeUncertaintyNanos.get(sv_index) * 0.299792458,2);
//            }else{
//                R[sv_index][sv_index] = Double.MAX_VALUE;
//            }
//        }
//        Matrix R_matrix = new Matrix(R);
//        Matrix R_matrix_inverse = R_matrix.inverse();

        int iterate_count = 0;
        double[][] user_position = new double[4][1];
        user_position[0][0] = roughPosition.ECEF_x;
        user_position[1][0] = roughPosition.ECEF_y;
        user_position[2][0] = roughPosition.ECEF_z;
        user_position[3][0] = roughPosition.clockBias;
        double[][] rou = new double[satelliteDataArrayList.size()][1];
        double[][] equation = new double[satelliteDataArrayList.size()][4];
        double[][] value = new double[satelliteDataArrayList.size()][1];
        Matrix delta_user_position_matrix,equation_matrix,value_matrix;
        Matrix user_position_matrix = new Matrix(user_position);
        do{
            iterate_count++;
            for(int sv_index = 0;sv_index < satelliteDataArrayList.size();++sv_index){
                rou[sv_index][0] = Math.sqrt(Math.pow(user_position[0][0]-satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_x,2) +
                        Math.pow(user_position[1][0]-satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_y,2) +
                        Math.pow(user_position[2][0]-satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_z,2));
                equation[sv_index][0] = (user_position[0][0] - satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_x) / rou[sv_index][0];
                equation[sv_index][1] = (user_position[1][0] - satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_y) / rou[sv_index][0];
                equation[sv_index][2] = (user_position[2][0] - satelliteDataArrayList.get(sv_index).satellitePosition.ECEF_z) / rou[sv_index][0];
                equation[sv_index][3] = 1;
                value[sv_index][0] = rou[sv_index][0] - satelliteDataArrayList.get(sv_index).getCorrectedPseudorange() + user_position[3][0];
            }
            equation_matrix = new Matrix(equation);
            value_matrix = new Matrix(value);
            //delta_user_position_matrix = equation_matrix.solve(value_matrix);
            delta_user_position_matrix = (equation_matrix.transpose().times(equation_matrix)).inverse().times(equation_matrix.transpose()).times(value_matrix);
            user_position_matrix = user_position_matrix.minus(delta_user_position_matrix);
            user_position = user_position_matrix.getArrayCopy();
        }while(delta_user_position_matrix.norm2() > 0.001 && iterate_count < 15);

        Position userPosition = new Position(user_position[0][0],user_position[1][0],user_position[2][0]);
        userPosition.setClockBias(user_position[3][0]);
        return userPosition;
    }


    /**
     *解开普勒方程(计算E_k用)
     */
    static public double keplerSolve(double e,double M,double tol){
        //double Mnorm = mod_2pi(M);
        double E0 = keplerInitial(e,M);
        double dE = tol + 1;
        int count = 0;
        double E = 0;
        while(dE > tol){
            E = E0 - keplerEps1(e,M,E0);
            dE = Math.abs(E - E0);
            E0 = E;
            ++count;
            if(count == 100){
                return -1;//易知解的范围在[0,2*pi],返回-1表示没有在100次循环内得到收敛的解
            }
        }
        return E;
    }
    static private double mod_2pi(double M){
        if(M < 0){
            while(M < 0){
                M += 2 * 3.1415927;
            }
        }else if(M > 2 * 3.1415927){
            while(M > 2 * 3.1415927) {
                M -= 2 * 3.1415927;
            }
        }
        return M;
    }
    static private double keplerInitial(double e,double M){
        double e_2 = e * e;
        double e_3 = e_2 * e;
        double cos_M = Math.cos(M);
        double E0 = M + (-0.5*e_3+e+(e_2+1.5*e_3*cos_M)*cos_M)*Math.sin(M);
        return E0;
    }
    static private double keplerEps1(double e,double M,double x){
        double eps1 = (x - e * Math.sin(x) - M)/(1 - e * Math.cos(x));
        return eps1;
    }


    static public double[] dateOfGPST(double GPST) {
        double secOfDay = GPST % (3600 * 24);
        double secOfHour = secOfDay % 3600;
        double sec = secOfHour % 60;
        double min = floor(secOfHour / 60);
        double hour = floor(secOfDay / 3600);
        double nDay = floor(GPST / (3600 * 24));
        double year = 1980;
        double month = 1;
        double day, yearDay;
        if (nDay >= 361) {
            day = 1;
            year = year + 1;
            nDay = nDay - 361;
            if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) {
                yearDay = 366;
            } else {
                yearDay = 365;
            }

            while (nDay >= yearDay) {
                nDay = nDay - yearDay;
                year = year + 1;
                if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) {
                    yearDay = 366;
                } else {
                    yearDay = 365;
                }
            }
        } else {
            day = 6;
        }

        double[] dayInMonth;
        if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) {
            dayInMonth = new double[]{31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        } else {
            dayInMonth = new double[]{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        }
        while (nDay >= dayInMonth[(int) month]) {
            nDay = nDay - dayInMonth[(int) month];
            month = month + 1;
        }
        day = day + nDay;

        return new double[]{year,month,day,hour,min,sec};

    }
}
