using UnityEngine;

public class Plane : MonoBehaviour
{
    public class Longitude
    {
        public float x, z, theta;
        public float U, W, Q;
        public float T;
        public float de;

        // Longitudinal derivatives
        public float Xu = -5.2172f;
        public float Xw = 122.39f;
        public float Zu = -280.43f;
        public float Zw = -794.08f;
        public float Zq = -198.81f;
        public float Mu = -3.2135f;
        public float Mw = -381.09f;
        public float Mq = -1413.2f;

        // Control derivatives 
        public float Xde = -3.5977f;
        public float Zde = 306.44f;
        public float Mde = 1647.6f;


    }
    public class Lateral
    {
        public float y, psi, phi;
        public float V, P, R;
        public float dr;

        // Lateral derivatives

        public float Yv = -55.02f;
        public float Yp = -1140.2f;
        public float Yr = 605.54f;
        public float Lv = -1171.8f;
        public float Lp = -54024f;
        public float Lr = 18642f;
        public float Nv = 21.28f;
        public float Np = -9323.9f;
        public float Nr = -1360.9f;
    }


    // public
    [System.NonSerialized] public float Airspeed = 0.000f; // Airspeed [m/s]
    [System.NonSerialized] public float alpha = 0.000f; // Angle of attack [deg]
    [System.NonSerialized] public float beta = 0.000f; // Side slip angle [deg]
    [System.NonSerialized] public float de = 0.000f; // Elevator angle [deg]
    [System.NonSerialized] public float dr = 0.000f; // Rudder angle [deg]
    [System.NonSerialized] public float dh = 0.000f; // Movement of c.g. [-]
    [System.NonSerialized] public float LocalGustMag = 0.000f; // Magnitude of local gust [m/s]
    [System.NonSerialized] public float LocalGustDirection = 0.000f; // Magnitude of local gust [m/s]
    [System.NonSerialized] public float nz = 0.000f; // Load factor [-]
    // Phisics
    static private float rho = 1.164f;
    static private float hE0 = 10.500f; // Altitude at Take-off [m]
    // At Cruise without Ground Effect
    static private float Airspeed0; // Magnitude of ground speed [m/s]
    static private float alpha0; // Angle of attack [deg]
    static private float CDp0; // Parasitic drag [-]
    static private float Cmw0; // Pitching momentum [-]
    static private float CLMAX; // Lift Coefficient [-]
    static private float CL0 = 0.000f; // Lift Coefficient [-]
    static private float CLw0 = 0.000f; // Lift Coefficient [-]
    static private float CLt0 = 0.000f; // Tail Coefficient [-]
    static private float epsilon0 = 0.000f; // Downwash
    // Plane
    static bool Downwash; // Conventional Tail: True, T-Tail: False
    static private float CL = 0.000f; // Lift Coefficient [-]
    static private float CD = 0.000f; // Drag Coefficient [-]
    static private float Cx = 0.000f; // X Force Coefficient [-]
    static private float Cy = 0.000f; // Y Force Coefficient [-]
    static private float Cz = 0.000f; // Z Force Coefficient [-]
    static private float Cl = 0.000f; // Rolling momentum [-]
    static private float Cm = 0.000f; // Pitching momentum [-]
    static private float Cn = 0.000f; // Yawing momentum [-]
    static private float dh0 = 0.000f; // Initial Mouse Position
    // Wing
    static private float Sw; // Wing area of wing [m^2]
    static private float bw; // Wing span [m]
    static private float cMAC; // Mean aerodynamic chord [m]
    static public float aw; // Wing Lift Slope [1/deg]
    static private float hw; // Length between Wing a.c. and c.g. [-]
    static private float AR; // Aspect Ratio [-]
    static private float ew; // Wing efficiency [-]
    static private float CLw = 0.000f; // Lift Coefficient [-]
    // Tail
    static private float St; // Wing area of tail [m^2]
    static private float at; // Tail Lift Slope [1/deg]
    static private float lt; // Length between Tail a.c. and c.g. [m]
    static private float VH; // Tail Volume [-]
    static private float deMAX; // Maximum elevator angle [deg]
    static private float tau; // Control surface angle of attack effectiveness [-]
    static private float CLt = 0.000f; // Lift Coefficient [-]
    // Fin
    static private float drMAX; // Maximum rudder angle
    // Ground Effect
    static private float CGEMIN; // Minimum Ground Effect Coefficient [-]
    static private float CGE = 0f; // Ground Effect Coefficient: CDiGE/CDi [-]
    // Stability derivatives
    static private float Cyb; // [1/deg]
    static private float Cyp; // [1/rad]
    static private float Cyr; // [1/rad]
    static private float Cydr; // [1/deg]
    static private float Cnb; // [1/deg]
    static private float Cnp; // [1/rad]
    static private float Cnr; // [1/rad]
    static private float Cndr; // [1/deg]
    static private float Clb; // [1/deg]
    static private float Clp; // [1/rad]
    static private float Clr; // [1/rad]
    static private float Cldr; // [1/deg]
    // Gust
    static private Vector3 Gust = Vector3.zero; // Gust [m/s]
    // Rotation
    static private float phi; // [deg]
    static private float theta;  // [deg]
    static private float psi; // [deg]

    public Lateral lateral = new Lateral();

    public Longitude longitude = new Longitude();

    private Rigidbody PlaneRigidbody;

    public float U0 = 12.0f;

    void Start()
    {
        // Get rigidbody component
        PlaneRigidbody = this.GetComponent<Rigidbody>();

        // Input Specifications
        InputSpecifications("QX-20");

        // Set take-off speed
        //if (GManager.instance.FlightMode == "BirdmanRally")
        //{
        //    //GManager.instance.Airspeed_TO = 5.0f; // Airspeed at take-off [m/s]
        //    PlaneRigidbody.velocity = Vector3.zero;
        //}
        //else if (GManager.instance.FlightMode == "TestFlight")
        //{ // 
        PlaneRigidbody.velocity = new Vector3(
                Airspeed0 * Mathf.Cos(Mathf.Deg2Rad * alpha0),
                -Airspeed0 * Mathf.Sin(Mathf.Deg2Rad * alpha0),
                0f
            );
        //}

        // Calculate CL at cluise
        CL0 = (PlaneRigidbody.mass * Physics.gravity.magnitude) / (0.5f * rho * Airspeed0 * Airspeed0 * Sw);
        CLt0 = (Cmw0 + CL0 * hw) / (VH + (St / Sw) * hw);
        CLw0 = CL0 - (St / Sw) * CLt0;
        if (Downwash) { epsilon0 = (CL0 / (Mathf.PI * ew * AR)) * Mathf.Rad2Deg; }

        dh0 = Screen.height / 2f; // Initial Mouse Position

        GManager.instance.plane = this;
        GManager.instance.Longitude = longitude;
        GManager.instance.Lateral = lateral;


    }
    private void FixedUpdate()
    {

        // Velocity and AngularVelocity
        float u = transform.InverseTransformDirection(PlaneRigidbody.velocity).x;
        float v = -transform.InverseTransformDirection(PlaneRigidbody.velocity).z;
        float w = -transform.InverseTransformDirection(PlaneRigidbody.velocity).y;
        float p = -transform.InverseTransformDirection(PlaneRigidbody.angularVelocity).x * Mathf.Rad2Deg;
        float q = transform.InverseTransformDirection(PlaneRigidbody.angularVelocity).z * Mathf.Rad2Deg;
        float r = transform.InverseTransformDirection(PlaneRigidbody.angularVelocity).y * Mathf.Rad2Deg;
        float hE = PlaneRigidbody.position.y+10.0f;

        // Force and Momentum
        Vector3 AerodynamicForce = Vector3.zero;
        Vector3 AerodynamicMomentum = Vector3.zero;
        Vector3 TakeoffForce = Vector3.zero;

        // Hoerner and Borst (Modified)
        CGE = (CGEMIN + 33f * Mathf.Pow((hE / bw), 1.5f)) / (1f + 33f * Mathf.Pow((hE / bw), 1.5f));

        // Get control surface angles
        //de = 0.000f;
        //dr = 0.000f;
        if (true)
        {
            dh = -(Input.mousePosition.y - dh0) * 0.0002f * 1.10f;
        }
        //Debug.Log(dh);

        longitude.de  = 0.0f * deMAX;
        lateral.dr = -0.0f * drMAX;

        if (Input.GetKey(KeyCode.LeftArrow))
        {
            lateral.dr += 10.0f;
        }else if (Input.GetKey(KeyCode.RightArrow))
        {
            lateral.dr -= 10.0f;
        }

        if (Input.GetKey(KeyCode.UpArrow))
        {
            longitude.de += 2.0f;
        }else if ((Input.GetKey(KeyCode.DownArrow)))
        {
            longitude.de += -2.0f;
        }

        if (Input.GetMouseButton(0)) { lateral.dr = drMAX; }
        else if (Input.GetMouseButton(1)) { lateral.dr = -drMAX; }

        //Debug.Log(dr);

        // Gust
        LocalGustMag = 1.0f * Mathf.Pow((hE / hE0), 1f / 7f);
        Gust = Quaternion.AngleAxis(0, Vector3.up) * (Vector3.right * LocalGustMag);
        Vector3 LocalGust = this.transform.InverseTransformDirection(Gust);
        float ug = 0 + 1e-10f;
        float vg = -0;
        float wg = -0;
        if (ug > 0) { LocalGustDirection = Mathf.Atan(vg / (ug + 1e-10f)) * Mathf.Rad2Deg; }
        else { LocalGustDirection = Mathf.Atan(vg / (ug + 1e-10f)) * Mathf.Rad2Deg + vg / Mathf.Abs((vg + 1e-10f)) * 180; }

        // Calculate angles
        Airspeed = Mathf.Sqrt((u + ug) * (u + ug) + (v + vg) * (v + vg) + (w + wg) * (w + wg));
        alpha = Mathf.Atan((w + wg) / (u + ug)) * Mathf.Rad2Deg;
        beta = Mathf.Atan((v + vg) / Airspeed) * Mathf.Rad2Deg;

        // Wing and Tail
        CLw = CLw0 + aw * (alpha - alpha0);
        CLt = CLt0 + at * ((alpha - alpha0) + (1f - CGE * (CLw / CLw0)) * epsilon0 + longitude.de * tau + ((lt - dh * cMAC) / Airspeed) * q);
        if (Mathf.Abs(CLw) > CLMAX) { CLw = (CLw / Mathf.Abs(CLw)) * CLMAX; } // Stall
        if (Mathf.Abs(CLt) > CLMAX) { CLt = (CLt / Mathf.Abs(CLt)) * CLMAX; } // Stall

        // Lift and Drag
        CL = CLw + (St / Sw) * CLt; // CL        
        CD = CDp0 * (1f + Mathf.Abs(Mathf.Pow((alpha / 9f), 3f))) + ((CL * CL) / (Mathf.PI * ew * AR)) * CGE; // CD

        // Force
        Cx = CL * Mathf.Sin(Mathf.Deg2Rad * alpha) - CD * Mathf.Cos(Mathf.Deg2Rad * alpha); // Cx       
        Cy = Cyb * beta + Cyp * (1f / Mathf.Rad2Deg) * ((p * bw) / (2f * Airspeed)) + Cyr * (1f / Mathf.Rad2Deg) * ((r * bw) / (2f * Airspeed)) + Cydr * lateral.dr; // Cy       
        Cz = -CL * Mathf.Cos(Mathf.Deg2Rad * alpha) - CD * Mathf.Sin(Mathf.Deg2Rad * alpha); // Cz


        // Torque
        Cl = Clb * beta + Clp * (1f / Mathf.Rad2Deg) * ((p * bw) / (2f * Airspeed)) + Clr * (1f / Mathf.Rad2Deg) * ((r * bw) / (2f * Airspeed)) + Cldr * lateral.dr; // Cl        
        Cm = Cmw0 + CLw * hw - VH * CLt + CL * dh; // Cm       
        Cn = Cnb * beta + Cnp * (1f / Mathf.Rad2Deg) * ((p * bw) / (2f * Airspeed)) + Cnr * (1f / Mathf.Rad2Deg) * ((r * bw) / (2f * Airspeed)) + Cndr * lateral.dr; // Cn

        Debug.Log(Cm);

        AerodynamicForce.x = 0.5f * rho * Airspeed * Airspeed * Sw * Cx;
        AerodynamicForce.y = 0.5f * rho * Airspeed * Airspeed * Sw * (-Cz);
        AerodynamicForce.z = 0.5f * rho * Airspeed * Airspeed * Sw * (-Cy);

        AerodynamicMomentum.x = 0.5f * rho * Airspeed * Airspeed * Sw * bw * (-Cl);
        AerodynamicMomentum.y = 0.5f * rho * Airspeed * Airspeed * Sw * bw * Cn;
        AerodynamicMomentum.z = 0.5f * rho * Airspeed * Airspeed * Sw * cMAC * Cm;


        float Distance = (PlaneRigidbody.position - GManager.instance.PlatformPosition).magnitude - 10f;
        if (GManager.instance.FlightMode == "BirdmanRally" && Distance < -0.5f)
        {

            CalculateRotation();

            float W = PlaneRigidbody.mass * Physics.gravity.magnitude;
            float L = 0.5f * rho * Airspeed * Airspeed * Sw * (Cx * Mathf.Sin(Mathf.Deg2Rad * theta) - Cz * Mathf.Cos(Mathf.Deg2Rad * theta));
            float N = (W - L) * Mathf.Cos(Mathf.Deg2Rad * 3.5f); // N=(W-L)*cos(3.5deg)
            float P = (PlaneRigidbody.mass * 0 * 0) / (2f * 10f); // P=m*Vto*Vto/2*L

            TakeoffForce.x = P;
            TakeoffForce.y = N * Mathf.Cos(Mathf.Deg2Rad * 3.5f);
            TakeoffForce.z = 0f;

            AerodynamicForce.z = 0f;
            AerodynamicMomentum.x = 0f;
            AerodynamicMomentum.y = 0f;
        }
        //Debug.Log(Distance);

        PlaneRigidbody.AddRelativeForce(AerodynamicForce, ForceMode.Force);
        PlaneRigidbody.AddRelativeTorque(AerodynamicMomentum, ForceMode.Force);
        PlaneRigidbody.AddForce(TakeoffForce, ForceMode.Force);

        nz = AerodynamicForce.y / (PlaneRigidbody.mass * Physics.gravity.magnitude);

        CalculateRotation();

        GManager.instance.de = de;
        GManager.instance.dr = dr;

        GManager.instance.Lateral.phi = phi;
        GManager.instance.Lateral.psi = psi;
        GManager.instance.Longitude.theta = theta;
    }


    void CalculateRotation()
    {
        float q1 =  transform.rotation.x;
        float q2 = -transform.rotation.y;
        float q3 = -transform.rotation.z;
        float q4 =  transform.rotation.w;
        float C11 = q1 * q1 - q2 * q2 - q3 * q3 + q4 * q4;
        float C22 = -q1 * q1 + q2 * q2 - q3 * q3 + q4 * q4;
        float C12 = 2f * (q1 * q2 + q3 * q4);
        float C13 = 2f * (q1 * q3 - q2 * q4);
        float C32 = 2f * (q2 * q3 - q1 * q4);

        phi = -Mathf.Atan(-C32 / C22) * Mathf.Rad2Deg;
        theta = -Mathf.Asin(C12) * Mathf.Rad2Deg;
        psi = -Mathf.Atan(-C13 / C11) * Mathf.Rad2Deg;
    }

    void InputSpecifications(string PlaneName)
    {
        if (PlaneName == "QX-20")
        {
            // Plane
            PlaneRigidbody.mass = 98.797f;
            PlaneRigidbody.centerOfMass = new Vector3(0f, 0.29f, 0f);
            PlaneRigidbody.inertiaTensor = new Vector3(1003f, 1045f, 58f);
            PlaneRigidbody.inertiaTensorRotation = Quaternion.AngleAxis(-9.112f, Vector3.forward);
            // Specification At Cruise without Ground Effect
            Airspeed0 = 9.600f; // Magnitude of ground speed [m/s]
            alpha0 = 1.459f; // Angle of attack [deg]
            CDp0 = 0.016f; // Parasitic drag [-]
            Cmw0 = -0.114f; // Pitching momentum [-]
            CLMAX = 1.700f;
            // Wing
            Sw = 18.816f; // Wing area of wing [m^2]
            bw = 26.679f; // Wing span [m]
            cMAC = 0.755f; // Mean aerodynamic chord [m]
            aw = 0.108f; // Wing Lift Slope [1/deg]
            hw = (0.323f - 0.250f); // Length between Wing a.c. and c.g.
            ew = 0.986f; // Wing efficiency
            AR = (bw * bw) / Sw; // Aspect Ratio
            // Tail
            Downwash = false; // Conventional Tail: True, T-Tail: False
            St = 1.526f; // Wing area of tail
            at = 0.088f; // Tail Lift Slope [1/deg]
            lt = 3.200f; // Length between Tail a.c. and c.g.
            deMAX = 10.000f; // Maximum elevator angle
            tau = 1.000f; // Control surface angle of attack effectiveness [-]
            VH = (St * lt) / (Sw * cMAC); // Tail Volume
            // Fin
            drMAX = 15.000f; // Maximum rudder angle            
            // Ground Effect
            CGEMIN = 0.293f; // Minimum Ground Effect Coefficient [-]
            // Stability derivatives
            Cyb = -0.003555f; // [1/deg]
            Cyp = -0.455493f; // [1/rad]
            Cyr = 0.143466f; // [1/rad]
            Cydr = 0.000888f; // [1/deg]
            Clb = -0.004049f; // [1/deg]
            Clp = -0.829690f; // [1/rad]
            Clr = 0.227736f; // [1/rad]
            Cldr = 0.000016f; // [1/deg]
            Cnb = -0.000500f; // [1/deg]
            Cnp = -0.132307f; // [1/rad]
            Cnr = 0.000942f; // [1/rad]
            Cndr = -0.000106f; // [1/deg]
        }
        else if (PlaneName == "ARG-2")
        {
            // Plane
            PlaneRigidbody.mass = 103.100f;
            PlaneRigidbody.centerOfMass = new Vector3(0f, -0.019f, 0f);
            PlaneRigidbody.inertiaTensor = new Vector3(961f, 1024f, 80f); //Ixx, Izz, Iyy
            PlaneRigidbody.inertiaTensorRotation = Quaternion.AngleAxis(-3.929f, Vector3.forward);
            // Specification At Cruise without Ground Effect
            Airspeed0 = 10.500f; // Magnitude of ground speed [m/s]
            alpha0 = 1.407f; // Angle of attack [deg] 
            CDp0 = 0.014f; // Parasitic drag [-]
            Cmw0 = -0.165f; // Pitching momentum [-]
            CLMAX = 1.700f;
            // Wing
            Sw = 18.009f; // Wing area of wing [m^2]
            bw = 23.350f; // Wing span [m]
            cMAC = 0.813f; // Mean aerodynamic chord [m]
            aw = 0.103f; // Wing Lift Slope [1/deg]
            hw = (0.3375f - 0.250f); // Length between Wing a.c. and c.g.
            ew = 0.986f; // Wing efficiency
            AR = (bw * bw) / Sw; // Aspect Ratio
            // Tail
            Downwash = true; // Conventional Tail: True, T-Tail: False
            St = 1.651f; // Wing area of tail
            at = 0.074f; // Tail Lift Slope [1/deg]
            lt = 3.200f; // Length between Tail a.c. and c.g.
            deMAX = 10.000f; // Maximum elevator angle
            tau = 1.000f; // Control surface angle of attack effectiveness [-]
            VH = (St * lt) / (Sw * cMAC); // Tail Volume
            // Fin
            drMAX = 15.000f; // Maximum rudder angle            
            // Ground Effect
            CGEMIN = 0.215f; // Minimum Ground Effect Coefficient [-]
            // Stability derivatives
            Cyb = -0.003764f; // [1/deg]
            Cyp = -0.411848f; // [1/rad]
            Cyr = 0.141631f; // [1/rad]
            Cydr = 0.001846f; // [1/deg]
            Clb = -0.003656f; // [1/deg]
            Clp = -0.816226f; // [1/rad]
            Clr = 0.219104f; // [1/rad]
            Cldr = 0.000032f; // [1/deg]
            Cnb = -0.000245f; // [1/deg]
            Cnp = -0.127263f; // [1/rad]
            Cnr = -0.002745f; // [1/rad]
            Cndr = -0.000308f; // [1/deg]
        }
        else if (PlaneName == "Meister23")
        {
            // Plane
            PlaneRigidbody.mass = 103.100f;
            PlaneRigidbody.centerOfMass = new Vector3(0f, -0.019f, 0f);
            PlaneRigidbody.inertiaTensor = new Vector3(961f, 1024f, 80f); //Ixx, Izz, Iyy
            PlaneRigidbody.inertiaTensorRotation = Quaternion.AngleAxis(-3.929f, Vector3.forward);
            // Specification At Cruise without Ground Effect
            Airspeed0 = 7.500f; // Magnitude of ground speed [m/s]
            alpha0 = 1.407f; // Angle of attack [deg] 
            CDp0 = 0.014f; // Parasitic drag [-]
            Cmw0 = -0.165f; // Pitching momentum [-]
            CLMAX = 1.700f;
            // Wing
            Sw = 18.009f; // Wing area of wing [m^2]
            bw = 23.350f; // Wing span [m]
            cMAC = 0.813f; // Mean aerodynamic chord [m]
            aw = 0.103f; // Wing Lift Slope [1/deg]
            hw = (0.3375f - 0.250f); // Length between Wing a.c. and c.g.
            ew = 0.986f; // Wing efficiency
            AR = (bw * bw) / Sw; // Aspect Ratio
            // Tail
            Downwash = true; // Conventional Tail: True, T-Tail: False
            St = 1.651f; // Wing area of tail
            at = 0.074f; // Tail Lift Slope [1/deg]
            lt = 3.200f; // Length between Tail a.c. and c.g.
            deMAX = 10.000f; // Maximum elevator angle
            tau = 1.000f; // Control surface angle of attack effectiveness [-]
            VH = (St * lt) / (Sw * cMAC); // Tail Volume
            // Fin
            drMAX = 15.000f; // Maximum rudder angle            
            // Ground Effect
            CGEMIN = 0.215f; // Minimum Ground Effect Coefficient [-]
            // Stability derivatives
            Cyb = -0.003450053f; // [1/deg]
            Cyp = -0.586670267f; // [1/rad]
            Cyr = 0.099732268f; // [1/rad]
            Cydr = 0.001846f; // [1/deg]
            Clb = -0.005119664f; // [1/deg]
            Clp = -1.042361891f; // [1/rad]
            Clr = 0.151214166f; // [1/rad]
            Cldr = 0.000032f; // [1/deg]
            Cnb = -0.000480888f; // [1/deg]
            Cnp = -0.084891343f; // [1/rad]
            Cnr = -0.002117846f; // [1/rad]
            Cndr = -0.000308f; // [1/deg]
        }
    }

}
