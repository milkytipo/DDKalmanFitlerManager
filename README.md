#DDKalmanFilterManager
------

Here I provide a class set to realize the real-time GNSS pseudorange double difference algorithm (RTD) in Java (actually for Android).

* The pseudorange algorithm comes from the paper **《MEMS-based IMU Assisted Real Time Difference Using Raw Measurements for Smartphone》**, I added it into this repo. ( You can neglect the IMU fusion way, because that fusion way is simple and a little stupid hhh). All GNSS raw data you can use [RawMeaLogger](https://github.com/milkytipo/RawMeaLogger) app to obtain(Notice:this app can obtain the raw measurement data from Android phone and output the Rinex 3.02 format. I suggest you take out the GNSS raw data process method in this app and leave other redundant methods out because there exist many chaotic methods I wrote.

* This class omitted many check items and for the detail procedure you can check in the [RTD_phone](https://github.com/milkytipo/RTD_phone). In this RTD repo, Initiallization is in the `main_PVT/ConfigReceiver_PVT/pvtCalculatorInitializing` and the Kalman Prediction is in the `main_PVT/pvt_forecast_filt/pvtEKF_init` and `main_PVT/pvt_forecast_filt/pvtEKF_prediction` and the Kalman Update procedure is in the `main_PVT/pointPos_LOG/kalmanPVT_resiraim_GPS/leastSquarePos_GPS_DD` and `main_PVT/pointPos_LOG/kalmanPVT_resiraim_GPS/kalmanPos_GPS1`. I suggest you can only read these .m files I mentioned above, because the RTD repo is too complicated, even myself dont want to see it again.
* If you want to use this class into your own Android App, I suggest you add the timestamp flag into every procedure in order to keep the sequence of every loop.
* I add two `.java` files, one you can directly copy it into your app, another is just for test to help you understand the functions of the class.(The test one contains the main function,in the main function, I initial the Update method's parameter value and you can just `javac` this .java and see the output of bur(the correct result should be bur[0] =140.76037,bur\[1]=145.643863,bur\[3]=-113.34793 ))


------
## Function and parameter illustration

Their are three primary method
1.`DDKalmanFilterManager()`
Constructor, initialize the P and Qw.
2. `Update(int[] activeChannel,double[][] satpos,double[][] satpos_ref,double[] obs,double[] obs_ref, PvtCalculator pvtCalculator,  double[] doppSmooth, double[] doppSmooth_ref, double[] elfore,double[] azfore,double[]cn0)`

`activeChannel   = new int[nmbOfSatellites]` // activeChannel array saves the prn of every GPS Satellite,nmbOfSatellites is the number of current valid satellites.
`satpos  =new double[6][32]` // 6 rows reserve the sat position_xyz and velocity_xyz; 32 volumns reserve all the 32 prn;`satpos_ref = new double[6][32] `
`obs = new double[32]`//every 32 GPS satellite pseudorange, invalid satellite ele is zer
`double[] obs_ref = new double[32];`
`doppSmooth =new double[32]` // //every 32 GPS satellite doppler value, invalid satellite ele is zer
`doppSmooth_ref =new double[32]`
`elfore = new double[32]` //every 32 GPS satellite elevation.invalid satellite ele is zero
`azfore = new double[32]` //every 32 GPS satellite elevation.invalid satellite ele is zero
`cn0 = new double[32]`//every 32 GPS satellite elevation.invalid satellite ele is zero
`PvtCalculator pvtCalculator`// memeber parameters are the Kalman state values.
public  class  PvtCalculator{
	public double[] stt_x  = new double[2];  //stt_*[0] refers the bur,namely baseline's length, stt_*[1] refers the velocity_xyz of user.
	public double[] stt_y = new double[2];
	public double[] stt_z  = new double[2];
	public double P[][] = new double[6][6];  //saves the convirance of every state.
};

3. `Predict(PvtCalculator pvc)`
pvc is a declared object,contains position xyz and velocity xyz and P&Qw in Kalman.
4. `getBur()` or `getBur(PvtCalculator pvc)`
return the bur xyz.

