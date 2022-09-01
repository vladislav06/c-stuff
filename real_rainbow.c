/*
                Colour Rendering of Spectra

                       by John Walker
                  http://www.fourmilab.ch/

                 Last updated: March 9, 2003

           This program is in the public domain.

    For complete information about the techniques employed in
    this program, see the World-Wide Web document:

             http://www.fourmilab.ch/documents/specrend/

    The xyz_to_rgb() function, which was wrong in the original
    version of this program, was corrected by:

            Andrew J. S. Hamilton 21 May 1999
            Andrew.Hamilton@Colorado.EDU
            http://casa.colorado.edu/~ajsh/

    who also added the gamma correction facilities and
    modified constrain_rgb() to work by desaturating the
    colour by adding white.

    A program which uses these functions to plot CIE
    "tongue" diagrams called "ppmcie" is included in
    the Netpbm graphics toolkit:
        http://netpbm.sourceforge.net/
    (The program was called cietoppm in earlier
    versions of Netpbm.)

*/

#include <math.h>
#include <stdio.h>

/* A colour system is defined by the CIE x and y coordinates of
   its three primary illuminants and the x and y coordinates of
   the white point. */

struct colourSystem {
  char *name;         /* Colour system name */
  double xRed, yRed,  /* Red x, y */
      xGreen, yGreen, /* Green x, y */
      xBlue, yBlue,   /* Blue x, y */
      xWhite, yWhite, /* White point x, y */
      gamma;          /* Gamma correction for system */
};

/* White point chromaticities. */

#define IlluminantC 0.3101, 0.3162         /* For NTSC television */
#define IlluminantD65 0.3127, 0.3291       /* For EBU and SMPTE */
#define IlluminantE 0.33333333, 0.33333333 /* CIE equal-energy illuminant */

/*  Gamma of nonlinear correction.

    See Charles Poynton's ColorFAQ Item 45 and GammaFAQ Item 6 at:

       http://www.poynton.com/ColorFAQ.html
       http://www.poynton.com/GammaFAQ.html

*/

#define GAMMA_REC709 0 /* Rec. 709 */

static struct colourSystem
    /* Name                  xRed    yRed    xGreen  yGreen  xBlue  yBlue White
       point        Gamma   */
    NTSCsystem = {"NTSC", 0.67, 0.33,        0.21,        0.71,
                  0.14,   0.08, IlluminantC, GAMMA_REC709},
    EBUsystem = {"EBU (PAL/SECAM)", 0.64,        0.33, 0.29, 0.60, 0.15, 0.06,
                 IlluminantD65,     GAMMA_REC709},
    SMPTEsystem = {"SMPTE", 0.630, 0.340,         0.310,       0.595,
                   0.155,   0.070, IlluminantD65, GAMMA_REC709},
    HDTVsystem = {"HDTV", 0.670, 0.330,         0.210,       0.710,
                  0.150,  0.060, IlluminantD65, GAMMA_REC709},
    CIEsystem = {"CIE",  0.7355, 0.2645,      0.2658,      0.7243,
                 0.1669, 0.0085, IlluminantE, GAMMA_REC709},
    Rec709system = {"CIE REC 709", 0.64, 0.33,          0.30,        0.60,
                    0.15,          0.06, IlluminantD65, GAMMA_REC709};

/*                          UPVP_TO_XY

    Given 1976 coordinates u', v', determine 1931 chromaticities x, y

*/

void upvp_to_xy(double up, double vp, double *xc, double *yc) {
  *xc = (9 * up) / ((6 * up) - (16 * vp) + 12);
  *yc = (4 * vp) / ((6 * up) - (16 * vp) + 12);
}

/*                          XY_TO_UPVP

    Given 1931 chromaticities x, y, determine 1976 coordinates u', v'

*/

void xy_to_upvp(double xc, double yc, double *up, double *vp) {
  *up = (4 * xc) / ((-2 * xc) + (12 * yc) + 3);
  *vp = (9 * yc) / ((-2 * xc) + (12 * yc) + 3);
}

/*                             XYZ_TO_RGB

    Given an additive tricolour system CS, defined by the CIE x
    and y chromaticities of its three primaries (z is derived
    trivially as 1-(x+y)), and a desired chromaticity (XC, YC,
    ZC) in CIE space, determine the contribution of each
    primary in a linear combination which sums to the desired
    chromaticity.  If the  requested chromaticity falls outside
    the Maxwell  triangle (colour gamut) formed by the three
    primaries, one of the r, g, or b weights will be negative.

    Caller can use constrain_rgb() to desaturate an
    outside-gamut colour to the closest representation within
    the available gamut and/or norm_rgb to normalise the RGB
    components so the largest nonzero component has value 1.

*/

void xyz_to_rgb(struct colourSystem *cs, double xc, double yc, double zc,
                double *r, double *g, double *b) {
  double xr, yr, zr, xg, yg, zg, xb, yb, zb;
  double xw, yw, zw;
  double rx, ry, rz, gx, gy, gz, bx, by, bz;
  double rw, gw, bw;

  xr = cs->xRed;
  yr = cs->yRed;
  zr = 1 - (xr + yr);
  xg = cs->xGreen;
  yg = cs->yGreen;
  zg = 1 - (xg + yg);
  xb = cs->xBlue;
  yb = cs->yBlue;
  zb = 1 - (xb + yb);

  xw = cs->xWhite;
  yw = cs->yWhite;
  zw = 1 - (xw + yw);

  /* xyz -> rgb matrix, before scaling to white. */

  rx = (yg * zb) - (yb * zg);
  ry = (xb * zg) - (xg * zb);
  rz = (xg * yb) - (xb * yg);
  gx = (yb * zr) - (yr * zb);
  gy = (xr * zb) - (xb * zr);
  gz = (xb * yr) - (xr * yb);
  bx = (yr * zg) - (yg * zr);
  by = (xg * zr) - (xr * zg);
  bz = (xr * yg) - (xg * yr);

  /* White scaling factors.
     Dividing by yw scales the white luminance to unity, as conventional. */

  rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
  gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
  bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;

  /* xyz -> rgb matrix, correctly scaled to white. */

  rx = rx / rw;
  ry = ry / rw;
  rz = rz / rw;
  gx = gx / gw;
  gy = gy / gw;
  gz = gz / gw;
  bx = bx / bw;
  by = by / bw;
  bz = bz / bw;

  /* rgb of the desired point */

  *r = (rx * xc) + (ry * yc) + (rz * zc);
  *g = (gx * xc) + (gy * yc) + (gz * zc);
  *b = (bx * xc) + (by * yc) + (bz * zc);
}

/*                            INSIDE_GAMUT

     Test whether a requested colour is within the gamut
     achievable with the primaries of the current colour
     system.  This amounts simply to testing whether all the
     primary weights are non-negative. */

int inside_gamut(double r, double g, double b) {
  return (r >= 0) && (g >= 0) && (b >= 0);
}

/*                          CONSTRAIN_RGB

    If the requested RGB shade contains a negative weight for
    one of the primaries, it lies outside the colour gamut
    accessible from the given triple of primaries.  Desaturate
    it by adding white, equal quantities of R, G, and B, enough
    to make RGB all positive.  The function returns 1 if the
    components were modified, zero otherwise.

*/

int constrain_rgb(double *r, double *g, double *b) {
  double w;

  /* Amount of white needed is w = - min(0, *r, *g, *b) */

  w = (0 < *r) ? 0 : *r;
  w = (w < *g) ? w : *g;
  w = (w < *b) ? w : *b;
  w = -w;

  /* Add just enough white to make r, g, b all positive. */

  if (w > 0) {
    *r += w;
    *g += w;
    *b += w;
    return 1; /* Colour modified to fit RGB gamut */
  }

  return 0; /* Colour within RGB gamut */
}

/*                          GAMMA_CORRECT_RGB

    Transform linear RGB values to nonlinear RGB values. Rec.
    709 is ITU-R Recommendation BT. 709 (1990) ``Basic
    Parameter Values for the HDTV Standard for the Studio and
    for International Programme Exchange'', formerly CCIR Rec.
    709. For details see

       http://www.poynton.com/ColorFAQ.html
       http://www.poynton.com/GammaFAQ.html
*/

void gamma_correct(const struct colourSystem *cs, double *c) {
  double gamma;

  gamma = cs->gamma;

  if (gamma == GAMMA_REC709) {
    /* Rec. 709 gamma correction. */
    double cc = 0.018;

    if (*c < cc) {
      *c *= ((1.099 * pow(cc, 0.45)) - 0.099) / cc;
    } else {
      *c = (1.099 * pow(*c, 0.45)) - 0.099;
    }
  } else {
    /* Nonlinear colour = (Linear colour)^(1/gamma) */
    *c = pow(*c, 1.0 / gamma);
  }
}

void gamma_correct_rgb(const struct colourSystem *cs, double *r, double *g,
                       double *b) {
  gamma_correct(cs, r);
  gamma_correct(cs, g);
  gamma_correct(cs, b);
}

/*                          NORM_RGB

    Normalise RGB components so the most intense (unless all
    are zero) has a value of 1.

*/

void norm_rgb(double *r, double *g, double *b) {
#define Max(a, b) (((a) > (b)) ? (a) : (b))
  double greatest = Max(*r, Max(*g, *b));

  if (greatest > 0) {
    *r /= greatest;
    *g /= greatest;
    *b /= greatest;
  }
#undef Max
}

/*                          SPECTRUM_TO_XYZ

    Calculate the CIE X, Y, and Z coordinates corresponding to
    a light source with spectral distribution given by  the
    function SPEC_INTENS, which is called with a series of
    wavelengths between 380 and 780 nm (the argument is
    expressed in meters), which returns emittance at  that
    wavelength in arbitrary units.  The chromaticity
    coordinates of the spectrum are returned in the x, y, and z
    arguments which respect the identity:

            x + y + z = 1.
*/

void wavelength_to_xyz(double wavelength, double *x, double *y, double *z) {
  int i;
  double lambda, X = 0, Y = 0, Z = 0, XYZ;

  /* CIE colour matching functions xBar, yBar, and zBar for
     wavelengths from 380 through 780 nanometers, every 5
     nanometers.  For a wavelength lambda in this range:

          cie_colour_match[(lambda - 380) / 5][0] = xBar
          cie_colour_match[(lambda - 380) / 5][1] = yBar
          cie_colour_match[(lambda - 380) / 5][2] = zBar

      To save memory, this table can be declared as floats
      rather than doubles; (IEEE) float has enough
      significant bits to represent the values. It's declared
      as a double here to avoid warnings about "conversion
      between floating-point types" from certain persnickety
      compilers. */

  static double cie_colour_match[][3] = {
      {0.0014, 0.0000, 0.0065}, {0.0015, 0.0000, 0.0070},
      {0.0016, 0.0000, 0.0077}, {0.0018, 0.0001, 0.0085},
      {0.0020, 0.0001, 0.0094}, {0.0022, 0.0001, 0.0105},
      {0.0025, 0.0001, 0.0120}, {0.0029, 0.0001, 0.0136},
      {0.0033, 0.0001, 0.0155}, {0.0037, 0.0001, 0.0177},
      {0.0042, 0.0001, 0.0201}, {0.0048, 0.0001, 0.0225},
      {0.0053, 0.0002, 0.0252}, {0.0060, 0.0002, 0.0284},
      {0.0068, 0.0002, 0.0320}, {0.0077, 0.0002, 0.0362},
      {0.0088, 0.0002, 0.0415}, {0.0100, 0.0003, 0.0473},
      {0.0113, 0.0003, 0.0536}, {0.0128, 0.0004, 0.0605},
      {0.0143, 0.0004, 0.0679}, {0.0156, 0.0004, 0.0741},
      {0.0171, 0.0005, 0.0810}, {0.0188, 0.0005, 0.0891},
      {0.0208, 0.0006, 0.0988}, {0.0232, 0.0006, 0.1102},
      {0.0263, 0.0007, 0.1249}, {0.0298, 0.0008, 0.1418},
      {0.0339, 0.0009, 0.1612}, {0.0384, 0.0011, 0.1830},
      {0.0435, 0.0012, 0.2074}, {0.0489, 0.0014, 0.2334},
      {0.0550, 0.0015, 0.2625}, {0.0618, 0.0017, 0.2949},
      {0.0693, 0.0019, 0.3311}, {0.0776, 0.0022, 0.3713},
      {0.0871, 0.0025, 0.4170}, {0.0976, 0.0028, 0.4673},
      {0.1089, 0.0031, 0.5221}, {0.1212, 0.0035, 0.5815},
      {0.1344, 0.0040, 0.6456}, {0.1497, 0.0046, 0.7201},
      {0.1657, 0.0052, 0.7980}, {0.1820, 0.0058, 0.8780},
      {0.1985, 0.0065, 0.9588}, {0.2148, 0.0073, 1.0391},
      {0.2299, 0.0081, 1.1141}, {0.2445, 0.0089, 1.1868},
      {0.2584, 0.0098, 1.2566}, {0.2716, 0.0107, 1.3230},
      {0.2839, 0.0116, 1.3856}, {0.2948, 0.0126, 1.4419},
      {0.3047, 0.0136, 1.4939}, {0.3136, 0.0146, 1.5414},
      {0.3216, 0.0157, 1.5844}, {0.3285, 0.0168, 1.6230},
      {0.3343, 0.0180, 1.6561}, {0.3391, 0.0192, 1.6848},
      {0.3430, 0.0204, 1.7094}, {0.3461, 0.0217, 1.7301},
      {0.3483, 0.0230, 1.7471}, {0.3496, 0.0243, 1.7599},
      {0.3501, 0.0256, 1.7695}, {0.3500, 0.0270, 1.7763},
      {0.3493, 0.0284, 1.7805}, {0.3481, 0.0298, 1.7826},
      {0.3464, 0.0313, 1.7833}, {0.3444, 0.0329, 1.7823},
      {0.3420, 0.0345, 1.7800}, {0.3392, 0.0362, 1.7765},
      {0.3362, 0.0380, 1.7721}, {0.3333, 0.0398, 1.7688},
      {0.3301, 0.0418, 1.7647}, {0.3267, 0.0438, 1.7593},
      {0.3229, 0.0458, 1.7525}, {0.3187, 0.0480, 1.7441},
      {0.3140, 0.0502, 1.7335}, {0.3089, 0.0526, 1.7208},
      {0.3033, 0.0550, 1.7060}, {0.2973, 0.0574, 1.6889},
      {0.2908, 0.0600, 1.6692}, {0.2839, 0.0626, 1.6473},
      {0.2766, 0.0653, 1.6226}, {0.2687, 0.0680, 1.5946},
      {0.2602, 0.0709, 1.5632}, {0.2511, 0.0739, 1.5281},
      {0.2406, 0.0770, 1.4849}, {0.2297, 0.0803, 1.4386},
      {0.2184, 0.0837, 1.3897}, {0.2069, 0.0872, 1.3392},
      {0.1954, 0.0910, 1.2876}, {0.1844, 0.0949, 1.2382},
      {0.1735, 0.0991, 1.1887}, {0.1628, 0.1034, 1.1394},
      {0.1523, 0.1079, 1.0904}, {0.1421, 0.1126, 1.0419},
      {0.1322, 0.1175, 0.9943}, {0.1226, 0.1226, 0.9474},
      {0.1133, 0.1279, 0.9015}, {0.1043, 0.1334, 0.8567},
      {0.0956, 0.1390, 0.8130}, {0.0873, 0.1446, 0.7706},
      {0.0793, 0.1504, 0.7296}, {0.0718, 0.1564, 0.6902},
      {0.0646, 0.1627, 0.6523}, {0.0580, 0.1693, 0.6162},
      {0.0519, 0.1763, 0.5825}, {0.0463, 0.1836, 0.5507},
      {0.0412, 0.1913, 0.5205}, {0.0364, 0.1994, 0.4920},
      {0.0320, 0.2080, 0.4652}, {0.0279, 0.2171, 0.4399},
      {0.0241, 0.2267, 0.4162}, {0.0207, 0.2368, 0.3939},
      {0.0175, 0.2474, 0.3730}, {0.0147, 0.2586, 0.3533},
      {0.0121, 0.2702, 0.3349}, {0.0099, 0.2824, 0.3176},
      {0.0079, 0.2952, 0.3014}, {0.0063, 0.3087, 0.2862},
      {0.0049, 0.3230, 0.2720}, {0.0037, 0.3385, 0.2588},
      {0.0029, 0.3548, 0.2464}, {0.0024, 0.3717, 0.2346},
      {0.0022, 0.3893, 0.2233}, {0.0024, 0.4073, 0.2123},
      {0.0029, 0.4256, 0.2010}, {0.0038, 0.4443, 0.1899},
      {0.0052, 0.4635, 0.1790}, {0.0070, 0.4830, 0.1685},
      {0.0093, 0.5030, 0.1582}, {0.0122, 0.5237, 0.1481},
      {0.0156, 0.5447, 0.1384}, {0.0195, 0.5658, 0.1290},
      {0.0240, 0.5870, 0.1201}, {0.0291, 0.6082, 0.1117},
      {0.0349, 0.6293, 0.1040}, {0.0412, 0.6502, 0.0968},
      {0.0480, 0.6707, 0.0901}, {0.0554, 0.6906, 0.0839},
      {0.0633, 0.7100, 0.0782}, {0.0716, 0.7280, 0.0733},
      {0.0805, 0.7453, 0.0687}, {0.0898, 0.7619, 0.0646},
      {0.0995, 0.7778, 0.0608}, {0.1096, 0.7932, 0.0573},
      {0.1202, 0.8082, 0.0539}, {0.1311, 0.8225, 0.0507},
      {0.1423, 0.8363, 0.0477}, {0.1538, 0.8495, 0.0449},
      {0.1655, 0.8620, 0.0422}, {0.1772, 0.8738, 0.0395},
      {0.1891, 0.8849, 0.0369}, {0.2011, 0.8955, 0.0344},
      {0.2133, 0.9054, 0.0321}, {0.2257, 0.9149, 0.0298},
      {0.2383, 0.9237, 0.0277}, {0.2511, 0.9321, 0.0257},
      {0.2640, 0.9399, 0.0238}, {0.2771, 0.9472, 0.0220},
      {0.2904, 0.9540, 0.0203}, {0.3039, 0.9602, 0.0187},
      {0.3176, 0.9660, 0.0172}, {0.3314, 0.9712, 0.0159},
      {0.3455, 0.9760, 0.0146}, {0.3597, 0.9803, 0.0134},
      {0.3741, 0.9841, 0.0123}, {0.3886, 0.9874, 0.0113},
      {0.4034, 0.9904, 0.0104}, {0.4183, 0.9929, 0.0095},
      {0.4334, 0.9950, 0.0087}, {0.4488, 0.9967, 0.0080},
      {0.4644, 0.9981, 0.0074}, {0.4801, 0.9992, 0.0068},
      {0.4960, 0.9998, 0.0062}, {0.5121, 1.0000, 0.0057},
      {0.5283, 0.9998, 0.0053}, {0.5447, 0.9993, 0.0049},
      {0.5612, 0.9983, 0.0045}, {0.5778, 0.9969, 0.0042},
      {0.5945, 0.9950, 0.0039}, {0.6112, 0.9926, 0.0036},
      {0.6280, 0.9897, 0.0034}, {0.6448, 0.9865, 0.0031},
      {0.6616, 0.9827, 0.0029}, {0.6784, 0.9786, 0.0027},
      {0.6953, 0.9741, 0.0026}, {0.7121, 0.9692, 0.0024},
      {0.7288, 0.9639, 0.0023}, {0.7455, 0.9581, 0.0022},
      {0.7621, 0.9520, 0.0021}, {0.7785, 0.9454, 0.0020},
      {0.7948, 0.9385, 0.0019}, {0.8109, 0.9312, 0.0019},
      {0.8268, 0.9235, 0.0018}, {0.8425, 0.9154, 0.0018},
      {0.8579, 0.9070, 0.0018}, {0.8731, 0.8983, 0.0017},
      {0.8879, 0.8892, 0.0017}, {0.9023, 0.8798, 0.0017},
      {0.9163, 0.8700, 0.0017}, {0.9298, 0.8598, 0.0016},
      {0.9428, 0.8494, 0.0016}, {0.9553, 0.8386, 0.0015},
      {0.9672, 0.8276, 0.0015}, {0.9786, 0.8163, 0.0014},
      {0.9894, 0.8048, 0.0013}, {0.9996, 0.7931, 0.0013},
      {1.0091, 0.7812, 0.0012}, {1.0181, 0.7692, 0.0012},
      {1.0263, 0.7570, 0.0011}, {1.0340, 0.7448, 0.0011},
      {1.0410, 0.7324, 0.0011}, {1.0471, 0.7200, 0.0010},
      {1.0524, 0.7075, 0.0010}, {1.0567, 0.6949, 0.0010},
      {1.0597, 0.6822, 0.0010}, {1.0617, 0.6695, 0.0009},
      {1.0628, 0.6567, 0.0009}, {1.0630, 0.6439, 0.0008},
      {1.0622, 0.6310, 0.0008}, {1.0608, 0.6182, 0.0008},
      {1.0585, 0.6053, 0.0007}, {1.0552, 0.5925, 0.0007},
      {1.0509, 0.5796, 0.0006}, {1.0456, 0.5668, 0.0006},
      {1.0389, 0.5540, 0.0005}, {1.0313, 0.5411, 0.0005},
      {1.0226, 0.5284, 0.0004}, {1.0131, 0.5157, 0.0004},
      {1.0026, 0.5030, 0.0003}, {0.9914, 0.4905, 0.0003},
      {0.9794, 0.4781, 0.0003}, {0.9665, 0.4657, 0.0003},
      {0.9529, 0.4534, 0.0003}, {0.9384, 0.4412, 0.0002},
      {0.9232, 0.4291, 0.0002}, {0.9072, 0.4170, 0.0002},
      {0.8904, 0.4050, 0.0002}, {0.8728, 0.3930, 0.0002},
      {0.8544, 0.3810, 0.0002}, {0.8349, 0.3689, 0.0002},
      {0.8148, 0.3568, 0.0002}, {0.7941, 0.3447, 0.0001},
      {0.7729, 0.3328, 0.0001}, {0.7514, 0.3210, 0.0001},
      {0.7296, 0.3094, 0.0001}, {0.7077, 0.2979, 0.0001},
      {0.6858, 0.2867, 0.0001}, {0.6640, 0.2757, 0.0001},
      {0.6424, 0.2650, 0.0000}, {0.6217, 0.2548, 0.0000},
      {0.6013, 0.2450, 0.0000}, {0.5812, 0.2354, 0.0000},
      {0.5614, 0.2261, 0.0000}, {0.5419, 0.2170, 0.0000},
      {0.5226, 0.2081, 0.0000}, {0.5035, 0.1995, 0.0000},
      {0.4847, 0.1911, 0.0000}, {0.4662, 0.1830, 0.0000},
      {0.4479, 0.1750, 0.0000}, {0.4298, 0.1672, 0.0000},
      {0.4121, 0.1596, 0.0000}, {0.3946, 0.1523, 0.0000},
      {0.3775, 0.1451, 0.0000}, {0.3608, 0.1382, 0.0000},
      {0.3445, 0.1315, 0.0000}, {0.3286, 0.1250, 0.0000},
      {0.3131, 0.1188, 0.0000}, {0.2980, 0.1128, 0.0000},
      {0.2835, 0.1070, 0.0000}, {0.2696, 0.1015, 0.0000},
      {0.2562, 0.0962, 0.0000}, {0.2432, 0.0911, 0.0000},
      {0.2307, 0.0863, 0.0000}, {0.2187, 0.0816, 0.0000},
      {0.2071, 0.0771, 0.0000}, {0.1959, 0.0728, 0.0000},
      {0.1852, 0.0687, 0.0000}, {0.1748, 0.0648, 0.0000},
      {0.1649, 0.0610, 0.0000}, {0.1554, 0.0574, 0.0000},
      {0.1462, 0.0539, 0.0000}, {0.1375, 0.0507, 0.0000},
      {0.1291, 0.0475, 0.0000}, {0.1212, 0.0446, 0.0000},
      {0.1136, 0.0418, 0.0000}, {0.1065, 0.0391, 0.0000},
      {0.0997, 0.0366, 0.0000}, {0.0934, 0.0342, 0.0000},
      {0.0874, 0.0320, 0.0000}, {0.0819, 0.0300, 0.0000},
      {0.0768, 0.0281, 0.0000}, {0.0721, 0.0263, 0.0000},
      {0.0677, 0.0247, 0.0000}, {0.0636, 0.0232, 0.0000},
      {0.0598, 0.0218, 0.0000}, {0.0563, 0.0205, 0.0000},
      {0.0529, 0.0193, 0.0000}, {0.0498, 0.0181, 0.0000},
      {0.0468, 0.0170, 0.0000}, {0.0437, 0.0159, 0.0000},
      {0.0408, 0.0148, 0.0000}, {0.0380, 0.0138, 0.0000},
      {0.0354, 0.0128, 0.0000}, {0.0329, 0.0119, 0.0000},
      {0.0306, 0.0111, 0.0000}, {0.0284, 0.0103, 0.0000},
      {0.0264, 0.0095, 0.0000}, {0.0245, 0.0088, 0.0000},
      {0.0227, 0.0082, 0.0000}, {0.0211, 0.0076, 0.0000},
      {0.0196, 0.0071, 0.0000}, {0.0182, 0.0066, 0.0000},
      {0.0170, 0.0061, 0.0000}, {0.0158, 0.0057, 0.0000},
      {0.0148, 0.0053, 0.0000}, {0.0138, 0.0050, 0.0000},
      {0.0129, 0.0047, 0.0000}, {0.0121, 0.0044, 0.0000},
      {0.0114, 0.0041, 0.0000}, {0.0106, 0.0038, 0.0000},
      {0.0099, 0.0036, 0.0000}, {0.0093, 0.0034, 0.0000},
      {0.0087, 0.0031, 0.0000}, {0.0081, 0.0029, 0.0000},
      {0.0076, 0.0027, 0.0000}, {0.0071, 0.0026, 0.0000},
      {0.0066, 0.0024, 0.0000}, {0.0062, 0.0022, 0.0000},
      {0.0058, 0.0021, 0.0000}, {0.0054, 0.0020, 0.0000},
      {0.0051, 0.0018, 0.0000}, {0.0047, 0.0017, 0.0000},
      {0.0044, 0.0016, 0.0000}, {0.0041, 0.0015, 0.0000},
      {0.0038, 0.0014, 0.0000}, {0.0036, 0.0013, 0.0000},
      {0.0033, 0.0012, 0.0000}, {0.0031, 0.0011, 0.0000},
      {0.0029, 0.0010, 0.0000}, {0.0027, 0.0010, 0.0000},
      {0.0025, 0.0009, 0.0000}, {0.0024, 0.0008, 0.0000},
      {0.0022, 0.0008, 0.0000}, {0.0020, 0.0007, 0.0000},
      {0.0019, 0.0007, 0.0000}, {0.0018, 0.0006, 0.0000},
      {0.0017, 0.0006, 0.0000}, {0.0015, 0.0006, 0.0000},
      {0.0014, 0.0005, 0.0000}, {0.0013, 0.0005, 0.0000},
      {0.0012, 0.0004, 0.0000}, {0.0012, 0.0004, 0.0000},
      {0.0011, 0.0004, 0.0000}, {0.0010, 0.0004, 0.0000},
      {0.0009, 0.0003, 0.0000}, {0.0009, 0.0003, 0.0000},
      {0.0008, 0.0003, 0.0000}, {0.0007, 0.0003, 0.0000},
      {0.0007, 0.0002, 0.0000}, {0.0006, 0.0002, 0.0000},
      {0.0006, 0.0002, 0.0000}, {0.0006, 0.0002, 0.0000},
      {0.0005, 0.0002, 0.0000}, {0.0005, 0.0002, 0.0000},
      {0.0004, 0.0002, 0.0000}, {0.0004, 0.0001, 0.0000},
      {0.0004, 0.0001, 0.0000}, {0.0004, 0.0001, 0.0000},
      {0.0003, 0.0001, 0.0000}, {0.0003, 0.0001, 0.0000},
      {0.0003, 0.0001, 0.0000}, {0.0003, 0.0001, 0.0000},
      {0.0003, 0.0001, 0.0000}, {0.0002, 0.0001, 0.0000},
      {0.0002, 0.0001, 0.0000}, {0.0002, 0.0001, 0.0000},
      {0.0002, 0.0001, 0.0000}, {0.0002, 0.0001, 0.0000},
      {0.0002, 0.0001, 0.0000}, {0.0002, 0.0001, 0.0000},
      {0.0001, 0.0001, 0.0000}, {0.0001, 0.0000, 0.0000},
      {0.0001, 0.0000, 0.0000}, {0.0001, 0.0000, 0.0000},
      {0.0001, 0.0000, 0.0000}, {0.0001, 0.0000, 0.0000},
      {0.0001, 0.0000, 0.0000}, {0.0001, 0.0000, 0.0000},
      {0.0001, 0.0000, 0.0000}, {0.0001, 0.0000, 0.0000},
      {0.0001, 0.0000, 0.0000}, {0.0001, 0.0000, 0.0000},
      {0.0001, 0.0000, 0.0000}, {0.0001, 0.0000, 0.0000},
      {0.0001, 0.0000, 0.0000}, {0.0001, 0.0000, 0.0000},
      {0.0000, 0.0000, 0.0000}, {0.0000, 0.0000, 0.0000},
      {0.0000, 0.0000, 0.0000}};

  X = cie_colour_match[(int)(wavelength - 380)][0];
  Y = cie_colour_match[(int)(wavelength - 380)][1];
  Z = cie_colour_match[(int)(wavelength - 380)][2];
  // printf("X:%f, Y:%f, Z:%f \n", X, Y, Z);
  XYZ = (X + Y + Z);
  *x = X / XYZ;
  *y = Y / XYZ;
  *z = Z / XYZ;
}

/*                            BB_SPECTRUM

    Calculate, by Planck's radiation law, the emittance of a black body
    of temperature bbTemp at the given wavelength (in metres).  */

double bbTemp = 5000; /* Hidden temperature argument
                         to BB_SPECTRUM. */
double bb_spectrum(double wavelength) {
  double wlm = wavelength * 1e-9; /* Wavelength in meters */

  return (3.74183e-16 * pow(wlm, -5.0)) /
         (exp(1.4388e-2 / (wlm * bbTemp)) - 1.0);
}
double white(double wavelength) { return 1; }

int main() {
  double t, x, y, z, r, g, b;
  struct colourSystem *cs = &SMPTEsystem;
  for (float lambda = 380; lambda < 750; lambda += 2) {
    wavelength_to_xyz(lambda, &x, &y, &z);
    xyz_to_rgb(cs, x, y, z, &r, &g, &b);
    constrain_rgb(&r, &g, &b);
    norm_rgb(&r, &g, &b);
    // printf("r:%.3f g:%.3f b:%.3f", r, g, b);
    printf("\033[48;2;%d;%d;%d m  \033[0m", (int)(r * 255), (int)(g * 255),
           (int)(b * 255));
    // printf("wavelength: %f\n",lambda);
  }

  return 0;
}
