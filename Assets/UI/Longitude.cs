using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Longtitude : MonoBehaviour
{
    Text text;
    // Start is called before the first frame update
    void Start()
    {
        text = GetComponent<Text>();
    }

    // Update is called once per frame
    void Update()
    {
        if (GManager.instance.Longitude != null)
        {
            text.text = string.Format("u={0} w={1}\nq={2}\nƒÆ={3}", GManager.instance.Longitude.U, GManager.instance.Longitude.W, GManager.instance.Longitude.Q, GManager.instance.Longitude.theta);
        }
    }
}
