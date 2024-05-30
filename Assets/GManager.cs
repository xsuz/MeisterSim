using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GManager : MonoBehaviour
{
    public static GManager instance;
    public double Distance { get; set; }
    public Vector3 Position { get; set; }

    public Vector3 PlatformPosition= Vector3.zero;

    public Plane.Longitude Longitude;

    public Plane.Lateral Lateral;

    public Plane plane;

    public float de = 0.0f;
    public float dr = 0.0f;

    public string FlightMode = "BirdmanRally";

    public GameObject Plane=null;
    // Start is called before the first frame update
    void Awake()
    {
        if(instance == null)
        {
            instance = this;
            Position = transform.position;
            DontDestroyOnLoad(this.gameObject);
        }
        else
        {
            Destroy(this.gameObject);
        }
    }
}
