using UnityEngine;
using System.Collections;

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
	//Movement control
        //public float speed, tilt;              // tilt factor
        //public Boundary boundary;       // movememnt boundary
		//public float smoothing = 5;     // this value is used for smoothing movement
		//private Vector3 smoothDirection;// used to smooth out mouse and touch control
		public Rigidbody myRigidbody;   // reference to rigitbody

        public GameObject shot;         // bullet prefab
        public Transform shotSpawn;     // the turret (bullet spawn location)
        public float fireRate = 0.5f;   
        private float nextFire = 0.0f;

		bool started = false;
		public bool Started { get { return started; } }

		//The current position and rotation of the player
		Vector3 r = new Vector3();//Unity’s default unit scale is 1 unit = 1 meter
		Quaternion q = new Quaternion {x=0, y=0, z=0, w=1};

		void NormalizeRotQ(ref Quaternion q) {
			float l2inv = 1.0f/Mathf.Sqrt (q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
			q.x *= l2inv;
			q.y *= l2inv;
			q.z *= l2inv;
			q.w *= l2inv;
		}

		public void VisualUpdate(Vector3 dr, Quaternion dtheta) {
			//Debug.Log ("visual update dr = " + dr);
			r += dr;
			q *= dtheta;
			NormalizeRotQ (ref q);
		}

		// Need a handle to the PiP camera transform, to unrotate it back to inertial frame
		public Camera pipCamera;
		Vector3 pipCamOffset;

        // Use this for initialization
        void Start()
        {
			Transform pipXform = pipCamera.transform;
			//Debug.Log ("PiP position = " + pipXform.position);
			pipCamOffset = pipXform.position;

            myRigidbody = GetComponent<Rigidbody>();//
			started = true;
        }

        void FixedUpdate()
        {
			#if KEYBOARD_MOUSE_INPUT
            // keyboard
            float moveHorizontal = Input.GetAxis("Horizontal");
            float moveVertical = Input.GetAxis("Vertical");

            // if keyboard direction key is pressed
            if (moveHorizontal != 0 || moveVertical != 0)
            {
                myRigidbody.velocity = new Vector3(moveHorizontal, 0.0f, moveVertical) * speed;

            }
            else
            {
                Vector3 pos = Input.mousePosition;

                pos.z = Camera.main.transform.position.y + 1;
                pos = Camera.main.ScreenToWorldPoint(pos);
                Vector3 origin = new Vector3(transform.position.x, transform.position.y, transform.position.z);//.zero;

                Vector2 currentPosition = new Vector3(pos.x, pos.z);
                Vector3 directionRaw = pos - origin;
                //Debug.Log("directionRaw.magnitude=" + directionRaw.magnitude);

                Vector3 direction = directionRaw.normalized;

                smoothDirection = Vector3.MoveTowards(smoothDirection, direction, smoothing);

                direction = smoothDirection;
                Vector3 movement = new Vector3(direction.x, 0, direction.z);
                myRigidbody.velocity = movement * speed;
                Debug.Log("currentPosition=" + currentPosition + "  Input.mousePosition=" + Input.mousePosition +
                    " pos =" + pos + " direction=" + direction + " smoothDirection=" + smoothDirection +
                    " movement=" + movement);

                //transform.position = Vector3.Lerp(transform.position, pos, Time.deltaTime* speed*2);
                //myRigidbody.velocity = new Vector3(moveHorizontal, 0.0f, moveVertical) * speed;                    
            }

            transform.position = new Vector3
            (
                Mathf.Clamp(transform.position.x, boundary.xMin, boundary.xMax),
                0.0f,
                Mathf.Clamp(transform.position.z, boundary.zMin, boundary.zMax)
            );
            myRigidbody.rotation = Quaternion.Euler(0, 0, myRigidbody.velocity.x * -tilt);
			#else
			//Debug.Log("Moving to " + this.r);
			myRigidbody.MovePosition(r); myRigidbody.MoveRotation(q);
			//myRigidbody.position = r; myRigidbody.rotation = q;
			pipCamera.transform.position = this.r + pipCamOffset;
			#endif
        }

		/*Don't fire a bullet in this navigation demo app
        void Update()
        {
            if ((Input.GetButton("Fire1") || Input.GetKeyDown(KeyCode.Space)) && Time.time > nextFire)
            {
                nextFire = Time.time + fireRate;
                Instantiate(shot, shotSpawn.position, shotSpawn.rotation);
                GetComponent<AudioSource>().Play();
            }
        }
        */
    }
}
