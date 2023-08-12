return 1;}


float3 ocean(float2 uv, float time, float steepness, float wavelength, float speed, float lacunarity, float2 direction,float iteration, float dir_rand,float noise_scale){
    const float PI = 3.14159265f;
    float3 offsetVector = float3(0,0,0);
    float3 tangent = float3(1,0,0);
    float3 binormal = float3(0,1,0);
    for(int i=0;i<iteration;i++){
        float k = 2.0f * PI / wavelength;
        int seed = wavelength/1000;
        int rotate_degrees = 360*noise((float2(i,i)+float2(0,0))*noise_scale)*float(i);
        //vec2 rotate(vec2 point, float degree, vec2 pivot)

        float2 random_dir = float2(1,0);
               random_dir = rotate(random_dir, rotate_degrees, float2(0,0));
               random_dir = normalize(random_dir)*dir_rand;
        float c = sqrt(9.8 / k);
        float2 d = normalize(direction);
       
        d = normalize(d+random_dir);
       
        float f = k * (dot(d,uv.xy) - c * time);
        
        float a = steepness / k;
        offsetVector.x += d.x * a * cos(f);
        offsetVector.z += a * sin(f);
        offsetVector.y += d.y * a * cos(f);

        tangent += (float3(
                    - d.x*d.x*steepness * sin(f),
                    - d.x*d.y*steepness * sin(f),
                      d.x    *steepness * cos(f)));
        binormal += float3(
                    - d.x * d.y * (steepness * sin(f)),
                    - d.y * d.y * (steepness * sin(f)),
                      d.y       * (steepness * cos(f)));
        wavelength /= lacunarity;
        steepness *=1.2;    
    }
    float3 normal = normalize(cross(binormal, tangent));

    return offsetVector;