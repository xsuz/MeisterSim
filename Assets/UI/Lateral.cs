using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Lateral : MonoBehaviour
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
        if (GManager.instance.Lateral != null )
            text.text = string.Format("v={0}\np={1} r={2}\nƒÕ={3}\nƒÓ={4}",
                GManager.instance.Lateral.V,
                GManager.instance.Lateral.P,
                GManager.instance.Lateral.R,
                GManager.instance.Lateral.psi,
                GManager.instance.Lateral.phi);
    }
}
