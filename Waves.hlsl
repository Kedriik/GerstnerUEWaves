return 1;}

float2 rotate(float2 p, float degree, float2 pivot) {
    float radAngle = -radians(degree);// "-" - clockwise
    float x = p.x;
    float y = p.y;

    float rX = pivot.x + (x - pivot.x) * cos(radAngle) - (y - pivot.y) * sin(radAngle);
    float rY = pivot.y + (x - pivot.x) * sin(radAngle) + (y - pivot.y) * cos(radAngle);

    return float2(rX, rY);
}

float3 mod289(float3 x) {
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

float2 mod289(float2 x) {
    return x - floor(x * (1.0 / 289.0)) * 289.0;
}

float3 permute(float3 x) {
    return mod289((x * 34.0 + 1.0) * x);
}

float3 taylorInvSqrt(float3 r) {
    return 1.79284291400159 - 0.85373472095314 * r;
}

// output noise is in range [-1, 1]
float snoise(float2 v) {
    const float4 C = float4(0.211324865405187,  // (3.0-sqrt(3.0))/6.0
                            0.366025403784439,  // 0.5*(sqrt(3.0)-1.0)
                            -0.577350269189626, // -1.0 + 2.0 * C.x
                            0.024390243902439); // 1.0 / 41.0

    // First corner
    float2 i  = floor(v + dot(v, C.yy));
    float2 x0 = v -   i + dot(i, C.xx);

    // Other corners
    float2 i1;
    i1.x = step(x0.y, x0.x);
    i1.y = 1.0 - i1.x;

    // x1 = x0 - i1  + 1.0 * C.xx;
    // x2 = x0 - 1.0 + 2.0 * C.xx;
    float2 x1 = x0 + C.xx - i1;
    float2 x2 = x0 + C.zz;

    // Permutations
    i = mod289(i); // Avoid truncation effects in permutation
    float3 p =
      permute(permute(i.y + float3(0.0, i1.y, 1.0))
                    + i.x + float3(0.0, i1.x, 1.0));

    float3 m = max(0.5 - float3(dot(x0, x0), dot(x1, x1), dot(x2, x2)), 0.0);
    m = m * m;
    m = m * m;

    // Gradients: 41 points uniformly over a line, mapped onto a diamond.
    // The ring size 17*17 = 289 is close to a multiple of 41 (41*7 = 287)
    float3 x = 2.0 * frac(p * C.www) - 1.0;
    float3 h = abs(x) - 0.5;
    float3 ox = floor(x + 0.5);
    float3 a0 = x - ox;

    // Normalise gradients implicitly by scaling m
    m *= taylorInvSqrt(a0 * a0 + h * h);

    // Compute final noise value at P
    float3 g = float3(
        a0.x * x0.x + h.x * x0.y,
        a0.y * x1.x + h.y * x1.y,
        g.z = a0.z * x2.x + h.z * x2.y
    );
    return 130.0 * dot(m, g);
}
float noise(float2 seed){
    return snoise(seed)*0.5 + 0.5;
}
float3 ocean_normal(float2 uv, float time, float steepness, float wavelength, float speed, float lacunarity, float2 direction,float iteration, float dir_rand, float noise_scale){
    const float PI = 3.14159265f;
    float3 offsetVector = float3(0,0,0);
    float3 tangent = float3(0,0,0);
    float3 binormal = float3(0,0,0);
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

        tangent += normalize(float3(
                    1 - d.x*d.x*steepness * sin(f),
                    - d.x*d.y*steepness*sin(f),
                    d.x*steepness  * cos(f)));
        binormal += float3(
                    -d.x * d.y * (steepness * sin(f)),
                    1 - d.y * d.y * (steepness * sin(f)),
                    d.y * (steepness * cos(f)));
        wavelength /= lacunarity;
        steepness *=1.2;    
    }
    float3 normal = normalize(cross(binormal, tangent));
    return -normal;