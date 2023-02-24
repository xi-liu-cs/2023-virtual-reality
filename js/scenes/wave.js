class complex
{
    constructor(x, y)
    {
        this.x = x;
        this.y = y;
    }

    check(a, b)
    {
        if(a == undefined || b == undefined) return false;
    }

    add(a, b)
    {
        if(!this.check(a, b)) return new complex(0, 0);
        return new complex(a.x + b.x, a.y + b.y);
    }

    sub(a, b)
    {
        if(!this.check(a, b)) return new complex(0, 0);
        return new complex(a.x - b.x, a.y - b.y);
    }

    mul(a, b)
    {
        if(!this.check(a, b)) return new complex(0, 0);
        return new complex(a.x * b.x - a.y * b.y, a.x * b.y + b.x * a.y);
    }

    scalar_mul(a, b)
    {
        if(!this.check(a, b)) return new complex(0, 0);
        return new complex(a.x * b, a.y * b);
    }

    scalar_div(a, b)
    {
        if(!this.check(a, b)) return new complex(0, 0);
        return new complex(a.x / b, a.y / b);
    }

    mod(a)
    {
        return Math.sqrt(a.x * a.x + a.y * a.y);
    }

    conj(a)
    {
        return new complex(a.x, -a.y);
    }

    polar(r, theta)
    {
        this.check(r, theta);
        return new complex(r * Math.cos(theta), r * Math.sin(theta));
    }
}

class fft
{
    constructor(n)
    {
        this.c = new complex();
        this.n = n;
    }

    fft(a)
    {
        let n = a.length,
        n_div_2 = n >> 1,
        out = new Array(n);
        if(n <= 1)
        {
            out[0] = a[0];
            if(a[0] && (Number.isNaN(a[0].x) || Number.isNaN(a[0].y)))
                out[0] = new complex(0.0, 0.0);
            return out;
        }
        let even = new Array(n_div_2),
        odd = new Array(n_div_2);
        for(let i = 0; i < n_div_2; ++i)
        {
            even[i] = a[2 * i];
            odd[i] = a[2 * i + 1];
        }
        even = this.fft(even);
        odd = this.fft(odd);
        let omega = -2.0 * Math.PI / n;
        for(let i = 0; i < n_div_2; ++i)
        {
            let polar = this.c.polar(1.0, i * omega);
            odd[i] = this.c.mul(odd[i], polar);
        }
        for(let i = 0; i < n_div_2; ++i)
        {
            out[i] = this.c.add(even[i], odd[i]);
            out[i + n_div_2] = this.c.sub(even[i], odd[i]);
        }
        return out;
    }
    
    inverse(a)
    {
        let n = a.length;
        for(let i = 0; i < n; ++i)
            a[i] = this.c.conj(a[i]);
        let out = this.fft(a);
        for(let i = 0; i < n; ++i)
            out[i] = this.c.conj(out[i]);
        return out;
    }

    ifft(a)
    {
        let n = this.n,
        n2 = n * n,
        a1 = new Array(n),
        a2 = new Array(n),
        a3 = new Array(n),
        out = new Array(n);
        for(let i = 0; i < n; ++i)
            a1[i] = this.inverse(a[i]);
        for(let i = 0; i < n; ++i)
        {
            a3[i] = new Array(n);
            for(let j = 0; j < n; ++j)
                a3[i][j] = this.c.scalar_div(a1[j][i], n2);
            a2[i] = this.inverse(a3[i]);
        }
        for(let i = 0; i < n; ++i)
        {
            out[i] = new Array(n);
            for (let j = 0; j < n; ++j)
                out[i][j] = a2[i][j].x;
        }
        return out; /* float out[n][n] */
    }
}

export function wave(u, v)
{
    let c = new complex();
    // Define constants
    const gravity = 9.81;
  
    // Define parameters
    const windSpeed = 10;
    const fetch = 200000;
  
    // Calculate Phillips spectrum
    const L = Math.pow(windSpeed, 2) / gravity;
    const w = 2 * Math.PI / L;
    const k = Math.sqrt(u * u + v * v) / w;
    const A = gravity / (w * w * Math.pow(k, 4));
    const B = Math.exp(-1 / (Math.pow(k * L, 2))) / Math.pow(k, 2);
    const spectrum = A * B;
  
    // Generate random complex numbers
    const n = 64;
    const r = new Array(n);
    for (let i = 0; i < n; i++) {
      r[i] = new complex(Math.random() - 0.5, Math.random() - 0.5);
    }
  
    // Perform fast Fourier transform
    const fftObj = new fft(n);
    const a = fftObj.fft(r);
  
    // Multiply Fourier coefficients by square root of spectrum
    for (let i = 0; i < n; i++) {
      const kx = (i - n / 2) * (2 * Math.PI / fetch);
      for (let j = 0; j < n; j++) {
        const ky = (j - n / 2) * (2 * Math.PI / fetch);
        const k = Math.sqrt(kx * kx + ky * ky);
        const index = i * n + j;
        if(a[index])
            a[index] = c.scalar_mul(a[index], Math.sqrt(spectrum / (k * k)));
      }
    }
  
    // Perform inverse Fourier transform
    const h_tilde = fftObj.ifft(a);
  
    // Calculate wave height
    let height = 0;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const index = i * n + j;
        height += h_tilde[i][j] * Math.exp(2 * Math.PI * complex.I * (u * (i - n / 2) + v * (j - n / 2)) / fetch);
      }
    }
  
    return [2 * u - 1, 2 * v - 1, height];
  }

  export function generate_ocean(n)
  {
    let L = 2000.0;
    let A = 20.0;
    let v = [];
    let k = [];
  
    for (let i = 0; i < n; i++) {
      v[i] = new Array(n).fill(new complex(0.0, 0.0));
      k[i] = new Array(n).fill(new complex(0.0, 0.0));
    }
  
    for (let i = 0; i < n; i++) {
      let kx = (i - n / 2) * (2 * Math.PI / L);
      for (let j = 0; j < n; j++) {
        let ky = (j - n / 2) * (2 * Math.PI / L);
        let mag = Math.sqrt(kx * kx + ky * ky);
        let index = i * n + j;
  
        if (mag < 0.000001) {
          v[i][j] = new complex(0.0, 0.0);
          k[i][j] = new complex(0.0, 0.0);
        } else {
          let phillips = A * Math.exp(-1.0 / (mag * L * mag * L)) / (mag * mag) * Math.pow(Math.abs(kx * Math.cos(45) + ky * Math.sin(45)), 2);
          let real = (Math.random() * 2 - 1) * Math.sqrt(phillips) / Math.sqrt(2);
          let imag = (Math.random() * 2 - 1) * Math.sqrt(phillips) / Math.sqrt(2);
  
          v[i][j] = new complex(real, imag);
          k[i][j] = new complex(kx, ky);
        }
      }
    }
  
    let v_fft = new fft(n);
    v = transpose(v);
    for (let i = 0; i < n; i++) {
      v[i] = v_fft.fft(v[i]);
      k[i] = v_fft.fft(k[i]);
      for (let j = 0; j < n; j++) {
        let len = v_fft.c.mod(k[i][j]);
        v[i][j] = v_fft.c.scalar_mul(v[i][j], 1.0 / Math.sqrt(len));
      }
    }
    v = transpose(v);
  
    let height = new Array(n);
    for (let i = 0; i < n; i++) {
      height[i] = new Array(n);
    }
  
    let v_ifft = new fft(n);
    v = transpose(v);
    for (let i = 0; i < n; i++) {
      v[i] = v_ifft.ifft(v[i]);
    }
    v = transpose(v);
  
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        height[i][j] = v[i][j].x;
      }
    }
  
    return height;
  }
  
  function transpose(matrix) {
    return matrix[0].map((col, i) => matrix.map(row => row[i]));
  }