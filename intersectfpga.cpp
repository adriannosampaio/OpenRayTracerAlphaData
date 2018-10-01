
/******************************
 * Extra information constants
 ******************************/


/**	Number of Vector-cores for triangle intersection calculation.
 * 	To change the number of triangle intersection Vector-cores
 * 	it's only necessary to update this constant
 */
#define NUM_CORES 8

/**	Helper constants
 **/
#define EPSILON 1.0e-2
#define INFINITY 1.0e+10


/** Triangle constants
 */
#define MAX_TRIS 50000
#define MAX_TDATA 450000
#define TRI_PARAMETERS 9

/* 	Ray constants
 *	Max resolution supported of 1920x1080 (Full HD)
 */
#define MAX_RAYS 2073600
replacing for 2000 rays
#define MAX_RDATA 12441600
#define RAY_PARAMETERS 6

/**	Helper macros that allow the use of define constants as arguments
 * 	for the HLS pragmas
 **/
#define PRAGMA_SUB(x) _Pragma (#x)
#define PRAGMA_HLS(x) PRAGMA_SUB(HLS x)

/**	Declare an HLS INTERFACE pragma directive that can receive define
 * 	constants as extra information
 *
 * 	Ex.: PRAGMA_INTERFACE(m_axi depth=MAX_TDATA port=i_tData)
 **/
#define PRAGMA_INTERFACE(x) PRAGMA_SUB(HLS INTERFACE x)
#define PRAGMA_AXI_LITE(PORT,BUNDLE) PRAGMA_SUB(HLS INTERFACE s_axilite port=PORT bundle=BUNDLE)
#define PRAGMA_AXI_FULL(PORT,DEPTH) PRAGMA_SUB(HLS INTERFACE m_axi depth=DEPTH port=PORT)

/**	Declare an HLS UNROLL pragma directive that can receive define
 * 	constants as extra information.
 *
 * 	Ex.: PRAGMA_UNROLL(factor=NUM_CORES)
 **/
//#define PRAGMA_UNROLL(FAC) PRAGMA_SUB(HLS UNROLL factor=FAC)
#define PRAGMA_UNROLL()    PRAGMA_SUB(HLS UNROLL)

/**	Declare an HLS ARRAY_PARTITION pragma directive that can receive define
 * 	constants as extra information.
 *
 * 	Ex.: PRAGMA_AP(variable=ids cyclic factor=NUM_CORES dim=1)
 **/

// pragma ARRAY_PARTITION of type "complete"
#define PRAGMA_AP_COMPLETE(VAR) \
	PRAGMA_SUB(HLS ARRAY_PARTITION variable=VAR complete dim=1)

// pragma ARRAY_PARTITION non-complete
#define PRAGMA_AP(VAR,TYPE,FACTOR) \
	PRAGMA_SUB(HLS ARRAY_PARTITION variable=VAR TYPE factor=FACTOR dim=1)


void intersectFPGA(
	int i_tNumber,           /// number of triangles
	volatile double *i_tData, /// triangle information.    Size = i_tNumber * Triangle attribute number(12)
	volatile int *i_tIds,    /// triangle id information. Size = i_tNumber
	int i_rNumber,           /// number of rays
	volatile double *i_rData, /// ray information          Size = i_rNumber * Ray Attr number(6)
	volatile int *o_tIds,     /// min intersection triangle ids for each ray. Size = i_rNumber
	volatile double *o_tIntersects   /// min intersection t value for each ray       Size = i_rNumber
	)
{

    // Allows to group the AXI-LITE function parameters in the
    // AXI-Lite control bundle
    PRAGMA_AXI_LITE(return, control)

    // AXI Pragmas for the triangle information input array
    PRAGMA_AXI_LITE(i_tData,control)
    PRAGMA_AXI_FULL(i_tData, MAX_TDATA)

    // AXI Pragmas for the triangle identifiers information input array
    PRAGMA_AXI_LITE(i_tIds,control)
    PRAGMA_AXI_FULL(i_tIds,MAX_TRIS)

    // AXI Pragmas for the ray information input array
    PRAGMA_AXI_LITE(i_rData,control)
    PRAGMA_AXI_FULL(i_rData,MAX_RDATA)

    // AXI Pragmas for the closest triangle id output array
    PRAGMA_AXI_LITE(o_tIds,control)
    PRAGMA_AXI_FULL(o_tIds,MAX_RAYS)

    // AXI Pragmas for the closest intersection point output array
    PRAGMA_AXI_LITE(o_tIntersects,control)
    PRAGMA_AXI_FULL(o_tIntersects,MAX_RAYS)

    // AXI Pragma for the number of triangles to be processed
    PRAGMA_AXI_LITE(i_tNumber,control)

    // AXI Pragma for the number of rays to be processed
    PRAGMA_AXI_LITE(i_rNumber,control)


    double ori[3], dir[3], p1[3], p2[3], p3[3];

    int ray, buff, core;
    rayLoop: for(ray = 0; ray < i_rNumber; ray++)
    {
   
        int rayBase = ray * TRI_PARAMETERS;

        /// getting ray data
        ori[0] = i_rData[rayBase],
        ori[1] = i_rData[rayBase + 1],
        ori[2] = i_rData[rayBase + 2];

        dir[0] = i_rData[rayBase + 3],
        dir[1] = i_rData[rayBase + 4],
        dir[2] = i_rData[rayBase + 5];

        double coreMinTriangleIntersections[NUM_CORES];
        PRAGMA_AP_COMPLETE(coreMinTriangleIntersections);

        int coreMinTriangleIds[NUM_CORES];
        PRAGMA_AP_COMPLETE(coreMinTriangleIds)

        inCoreMemset: for(core = 0; core < NUM_CORES; core++)
        {
            coreMinTriangleIntersections[core] = INFINITY;
            coreMinTriangleIds[core] = -1;
        }

        triBuffer: for(buff = 0; buff < i_tNumber; buff += NUM_CORES)
        {
            double tris[NUM_CORES * TRI_PARAMETERS];
            PRAGMA_AP(tris,block,NUM_CORES)

            int ids[NUM_CORES];
            PRAGMA_AP_COMPLETE(ids)

            memcpy(
                tris,
                i_tData + buff * TRI_PARAMETERS,
                TRI_PARAMETERS * NUM_CORES * sizeof(double)
            );

            memcpy(
                ids,
                i_tIds + buff,
                NUM_CORES * sizeof(int)
            );

            triLoop: for(core = 0; core < NUM_CORES; core++)
            {
                PRAGMA_UNROLL()

                int tri = buff + core;
                int triBase = core * TRI_PARAMETERS;

                p1[0] = tris[triBase],
                p1[1] = tris[triBase + 1],
                p1[2] = tris[triBase + 2];

                p2[0] = tris[triBase + 3],
                p2[1] = tris[triBase + 4],
                p2[2] = tris[triBase + 5];

                p3[0] = tris[triBase + 6],
                p3[1] = tris[triBase + 7],
                p3[2] = tris[triBase + 8];

                double a = p1[0] - p2[0], b = p1[0] - p3[0],
                   c = dir[0], d = p1[0] - ori[0];
                double e = p1[1] - p2[1], f = p1[1] - p3[1],
                   g = dir[1], h = p1[1] - ori[1];
                double i = p1[2] - p2[2], j = p1[2] - p3[2],
                   k = dir[2], l = p1[2] - ori[2];
                double m = f * k - g * j, n = h * k - g * l,
                   p = f * l - h * j;
                double q = g * i - e * k, s = e * j - f * i;
                double inv_denom = 1.0 / (a * m + b * q + c * s);
                double e1 = d*m - b*n - c*p;
                double beta = e1*inv_denom;

                if(beta >= 0.0)
                {
                    double r = e*l - h*i;
                    double e2 = a*n + d*q + c*r;
                    double gamma = e2 * inv_denom;

                    if(gamma >= 0.0)
                    {
                        if(beta + gamma <= 1.0)
                        {
                            double e3 = a*p - b*r + d*s;
                            double t = e3 * inv_denom;

                            if(t >= EPSILON)
                            {
                                if(t < coreMinTriangleIntersections[core]) /// if the intersection point is the closest yet
                                {
                                    /// sets the minimum t value of the ray
                                    coreMinTriangleIntersections[core] = t;
                                    
                                    /// and sets the id of the closest triangle
                                    coreMinTriangleIds[core] = ids[core];
                                }
                            }
                        }
                    } // gamma >= 0.0
                } // beta >= 0.0
            } // triLoop
        } // buffLoop

        double tmin = INFINITY;
        int idMin = -1;
        for(core = 0; core < NUM_CORES; core++)
        {
            if(coreMinTriangleIds[core] != -1 &&
            coreMinTriangleIntersections[core] < tmin) /// if no triangle was hit
            {
                /// sets output id of the ray to -1
                /// meaning that no triangle will be shaded
                tmin = coreMinTriangleIntersections[core];
                idMin = coreMinTriangleIds[core];
            }
        }
        
        o_tIds[ray] = idMin;
        o_tIntersects[ray] = tmin;

    }
}

