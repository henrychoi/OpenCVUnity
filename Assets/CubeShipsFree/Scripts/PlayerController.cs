using UnityEngine;
using System;
using UnityEngine.Assertions;
using UnityEngine.UI;

// Heavily based on Unity Space Shooter tutorial
namespace CubeSpaceFree
{

    [System.Serializable]
    public class Boundary
    {
        public float xMin, xMax, zMin, zMax;
    }

    public class PlayerController : MonoBehaviour
    {
		enum State {
			Resetting,
			Updating
		}
		State state = State.Updating;
		int nZVU;

		#region MOVEMENT CONTROL
        //public float speed, tilt;              // tilt factor
        //public Boundary boundary;       // movememnt boundary
		//public float smoothing = 5;     // this value is used for smoothing movement
		//private Vector3 smoothDirection;// used to smooth out mouse and touch control
		public Rigidbody myRigidbody;//Necessary for physics; will the player be subject to physics?
		#endregion

		public GameObject enemy;
		TextMesh enemy_text_mesh;

        public GameObject shot;         // bullet prefab
        Transform turret_xform;// the turret (bullet spawn location)
        public float fireRate = 0.5f;   
        // float nextFire = 0.0f;

		// Changed in Unity editor Edit -> Project Settings -> Time -> Fixed Timestep
		const int F_fixed_update = 50;//default is 50
		const float Ts = 1.0f / F_fixed_update;

		#region sensor stats
		Vector3 w_bias = new Vector3(); //gyro bias

		// accel specific gravity [g] window;
		//Specific acc INCLUDES gravity (i.e. REACTION to gravity NOT subtracted)
		Vector3Window a_inb_lb_ens = new Vector3Window ((int)(1.0f * F_fixed_update), null);//"acc");
		// gyro [dps] window
		Vector3Window w_inb_lb_ens = new Vector3Window ((int)(1.0f * F_fixed_update), null);
			//, m_ens = new Vector3Window (0.5 * F_fixed_update) // magnetic [uT] field window
		#endregion

		bool started = false;
		public bool Started { get { return started; } }

		public Text statusText;


		//The current position and rotation of the player in local tangtent frame (l)
		Vector3 r_inl_lb = new Vector3()
			, r_inU_lb  = new Vector3(); //same vector but scaled and rotated from physics to Unity
		Vector3Window r_inl_lb_ens = new Vector3Window (F_fixed_update, null); 

		static Quaternion Q_lb = new Quaternion (0, 0, Mathf.Sin (Mathf.PI / 4), Mathf.Cos (Mathf.PI / 4))
		                         * new Quaternion (Mathf.Sin (Mathf.PI / 4), 0, 0, Mathf.Cos (Mathf.PI / 4))
			//, Q2_lb = new Quaternion (0, Mathf.Sin (-Mathf.PI / 4), 0, Mathf.Cos (-Mathf.PI / 4))
			//	* new Quaternion (Mathf.Sin (Mathf.PI / 2), 0, 0, Mathf.Cos (Mathf.PI / 2))
			;
		Quaternion q_lb = Q_lb
			, q_Ub = Quaternion.identity
			, q_inc_lb = Quaternion.identity
		;

		static Quaternion Q_Ul = //Will be used often, so calculate once
			new Quaternion(0, Mathf.Sin(-Mathf.PI/4), 0, Mathf.Cos(-Mathf.PI/4))
			* new Quaternion(Mathf.Sin(-Mathf.PI/4), 0, 0, Mathf.Cos(-Mathf.PI/4))
			;

		public void VisualUpdate(Vector3 dr, Quaternion dtheta) {
			//Cutting the tie from visual odometry (not that optical flow is visual odometry)
//			r += dr * Scale_w2game;
//			q *= dtheta;
//			Rternion.Normalize (ref q);
		}

		// Need a handle to the PiP camera transform, to unrotate it back to inertial frame
		public Camera attitudeCam, positionCam;//These are connected to this class in the inspector
		Vector3 attitudeCamOffset, positionCamOffset;
		bool yIsUp = true;//alternative is up

		// Use this for initialization
        void Start()
		{
			//Debug.Log ("platform =" + Application.platform);
			if (Application.platform == RuntimePlatform.OSXEditor
			    || Application.platform == RuntimePlatform.WindowsEditor) {
				yIsUp = true;
				Debug.Log ("Using Unity Remote");
			}

			Input.compensateSensors = true;

			//Must be faster or equal rate than FixedUpdate interval
			Input.gyro.updateInterval = 1.0f/F_fixed_update;

			//Chase camera position offsets
			attitudeCamOffset = attitudeCam.transform.position;
			positionCamOffset = positionCam.transform.position;

            myRigidbody = GetComponent<Rigidbody>();
			Assert.IsNotNull (myRigidbody);

			turret_xform = transform.Find ("Turret");
			Assert.IsNotNull (turret_xform);

			enemy_text_mesh = enemy.GetComponent<TextMesh> ();
			Assert.IsNotNull (enemy_text_mesh);

			statusText.gameObject.SetActive(false);

			started = true;
        }


		//Vector3 origin = new Vector3();
		//void OnDrawGizmos() { //This IS called, but I just did not see the ray
		//	Debug.Log ("I am here");
		//	Debug.DrawRay(r, a, Color.yellow, 1.0f, false);
			//Ray ray = new Ray {}
			//Gizmos.DrawLine(origin, new Vector3 {x=1, y=0, z=0});
		//}
	    //const float HALF_PI = 0.5f * Mathf.PI;
	    //int case_num = 0;

		public void Reset()
		{
			//statusText.text = "Please hold still";
			statusText.gameObject.SetActive(true);
			nZVU = 0;
			state = State.Resetting;
		}
	
        void FixedUpdate()
        {   //Debug.Log("dT = " + Time.deltaTime);

			//not sure if this is compensated, but still estimate bias
			a_inb_lb_ens.Add(Input.acceleration);
			//Debug.Log(String.Format("acc = ({0:F},{1:F},{2:F})", acc.x, acc.y, acc.z));

			//I will do further bias estimation using the compensated gyro values
			w_inb_lb_ens.Add(Input.gyro.rotationRateUnbiased);

			//m_ens.Add(Input.compass.rawVector);
			//float heading = Input.compass.magneticHeading;//[deg]
			switch(state) {
			case State.Resetting:
				if (++nZVU < F_fixed_update)
					break;

				w_bias = w_inb_lb_ens.Average;
				//Coarse tip/tilt estimate, see Groves 2nd ed. section 5.6, p. 198
				Vector3 a_b = a_inb_lb_ens.Average.normalized; //specific force includes gravity

				//Calculate the nominal case
				float roll, pitch;
				if (yIsUp) {
					roll = Mathf.Asin (-a_b.x);
					pitch = Mathf.Atan2 (-a_b.z, -a_b.y);
				} else {
					roll = Mathf.Asin (a_b.y);
					pitch = Mathf.Atan2 (-a_b.z, -a_b.x);
				}
				//statusText.text = String.Format ("roll {0:N1} pitch {1:N1}, a_b {2}"
				//	, roll, pitch, a_b);

				//Assert.IsFalse(case_num == 2);
				//case_num = 1;

				//Reset quaternion to yaw degrees around the z axis,
				//I don't use Quaternion.Euler(Mathf.Rad2Deg * roll, Mathf.Rad2Deg * pitch, 0)
				//to avoid, thinking about left-handed rotation.
				roll *= 0.5f; // Halve the Euler angles feed to quaternion construction
				pitch *= 0.5f;

				Quaternion q = // Apply the roll first, then pitch
					new Quaternion (0, 0, Mathf.Sin (roll), Mathf.Cos (roll))
					* new Quaternion (Mathf.Sin (pitch), 0, 0, Mathf.Cos (pitch));
				//q_lb = Q_lb
				//	* Quaternion.Euler(pitch, 0, roll);//minor miracle: the order matches what I need
				//Debug.Log (String.Format ("roll {0}, pitch {1} a {2} q {3}", roll, pitch, a_b, q));

				const float THRESHOLD = 0.8f;
				if (Mathf.Abs (a_b.x) < THRESHOLD) { // 1st case: cos(roll) reasonable value
					q_lb = Q_lb * q;
					statusText.text = String.Format ("roll {0:N1} pitch {1:N1}, a_b {2}", roll, pitch, a_b);
				} else { // 2nd case, roll near +/- 90; ay << ax
					float yaw;
					if (yIsUp) {
						roll = -Mathf.Sign (a_b.x) * Mathf.Acos (-a_b.y);
					} else {
						roll = Mathf.Sign (a_b.y) * Mathf.Acos (-a_b.x);
					}
					yaw = Mathf.Asin (a_b.z / Mathf.Sin (roll));
					statusText.text = String.Format ("roll {0:N1} yaw {1:N1}, a_b {2}"
						, roll, yaw, a_b);

					//Assert.IsFalse(case_num == 1);
					//case_num = 2;

					//Reset quaternion to yaw degrees around the z axis,
					//I don't use Quaternion.Euler(Mathf.Rad2Deg * roll, Mathf.Rad2Deg * pitch, 0)
					//to avoid, thinking about left-handed rotation.
					roll *= 0.5f; // Halve the Euler angles feed to quaternion construction
					yaw *= 0.5f;

					Quaternion q2 =// Apply the roll first, then yaw about Y
						new Quaternion (0, 0, Mathf.Sin (roll), Mathf.Cos (roll))
						* new Quaternion (0, Mathf.Sin (yaw), 0, Mathf.Cos (yaw));
					//Debug.Log (String.Format("roll {0}, yaw {1} a {2} q {3}", roll, yaw, a_b, q_lb));

					q_lb = Q_lb * Quaternion.Slerp (q, q2, Mathf.Abs (a_b.x) - THRESHOLD);
				}
				//Coarse heading estimate, according to Groves 2nd ed. section 6.1.1, p. 219
				//			float s_roll = Mathf.Sin (roll), c_roll = Mathf.Cos (roll)
				//				, s_pitch = Mathf.Sin (pitch), c_pitch = Mathf.Cos (pitch)
				//			    , m_num = -m_b.y * c_roll + m_b.z * s_roll
				//				, m_den = m_b.x + m_b.y + m_b.z;

				r_inl_lb.Set (0, 0, 0);//Shall I reset the position at 1.5 m height?

				//statusText.gameObject.SetActive(false);
				//state = State.Updating;
				break;

			default: //navigation update
				//q_inc_lb.Set (0, 0, 0.5f * 0.01f, 1);//For now, let's make up a rotation vector
				//q_lb = q_lb * q_inc_lb;//Rotate the current attitude by the attitude increment
				//Note: Unity rotation sequence is left-to-right (the opposite of math operations)
				//Rternion.Normalize (ref q_lb);
				break;
			}

			// Update the Unity position and attitude
			//r_inU_lb.Set(-r_inl_lb.y, r_inl_lb.z, r_inl_lb.x);
			//const int Scale_w2game = 100;

			// physics is in 10 cm scale (phone dimension); Unity is in m
			r_inU_lb *= (1000/10);

			q_Ub = Q_Ul * q_lb; //Rotate the q_lb by Q_lU to go to Unity
			//q_Ub.Set(-q_Ub.x, -q_Ub.z, -q_Ub.y, q_Ub.w); //then go to Unity frame (LH)
			Rternion.Normalize(ref q_Ub);
			//Debug.Log("Moving to " + this.r);
			myRigidbody.MovePosition(r_inU_lb); myRigidbody.MoveRotation(q_Ub);
			//myRigidbody.position = r; myRigidbody.rotation = q;

			r_inl_lb_ens.Add(r_inl_lb);

			//Move the chase cameras (shown in PiP) with the player
			attitudeCam.transform.position = r_inl_lb + attitudeCamOffset;
			positionCam.transform.position = r_inl_lb_ens.Average + positionCamOffset;
        }

        void Update()
        {
			RaycastHit hitInfo;
			if (Physics.Raycast(r_inl_lb //The starting point of the ray in world coordinates.
				, turret_xform.forward // direction
				, out hitInfo //hitInfo will contain more information about where the collider was hit
				//, float maxDistance = Mathf.Infinity
				//, int layerMask = DefaultRaycastLayers
				//, QueryTriggerInteraction queryTriggerInteraction = QueryTriggerInteraction.UseGlobal
				)
			) { // can hit the enemy
				GameObject hit = hitInfo.transform.gameObject;
				if (hit == enemy) {
					enemy_text_mesh.color = Color.red;
					return;
				}
			}
			// else will not hit the target
			enemy_text_mesh.color = Color.black;

			/*Don't fire a bullet in this navigation demo app
            if ((Input.GetButton("Fire1") || Input.GetKeyDown(KeyCode.Space)) && Time.time > nextFire)
            {
                nextFire = Time.time + fireRate;
                Instantiate(shot, turret_xform.position, turret_xform.rotation);
                GetComponent<AudioSource>().Play();
            }
			        */
        }
    }
}
