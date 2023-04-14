/**
 * Some functions to help with 4x4 matrix manipulation
 * 
 * Most of the code is based on wikipedia and i tried to style the syntax like m4.js
 * Although the code is similar, no actual code was copied from m4.js
 * 
 */

var Neo4 = {
    // testMat: [2,0,-1,3,1,0,-2,0,7,2,16,24,-3,0,1,-6],

    /**
     *  Vec3 portion of grouping
     */
    vec3:{
        /**
         * constructs vec3 from values
         * @param {array[3]} vals 
         * @returns vec3
         */
        construct: function(vals){
            vals = vals || [];
            const v3 = new Float32Array(3);
            v3[0] = vals[0] || 0;
            v3[1] = vals[1] || 0;
            v3[2] = vals[2] || 0;
            return v3;
        },
        /**
         * Normalizes a given vector to 1
         * @param {vec3} v 
         * @param {vec3} dst 
         * @returns vec3
         */
        normalize: function(v, dst){
            dst = dst || new Float32Array(3);

            const v0 = v[0];
            const v1 = v[1];
            const v2 = v[2];
            const l = Math.sqrt(v0*v0 + v1*v1 + v2*v2);
            dst[0] = v0/l;
            dst[1] = v1/l;
            dst[2] = v2/l;

            return dst;
        },
        /**
         * adds two vec3s
         * @param {vec3} v0 
         * @param {vec3} v1 
         * @param {vec3} dst 
         * @returns vec3
         */
        add: function(v0, v1, dst){
            dst = dst || new Float32Array(3);

            dst[0] = v0[0] + v1[0];
            dst[1] = v0[1] + v1[1];
            dst[2] = v0[2] + v1[2];

            return dst;
        },

        /**
         * subtracts v1 from v0
         * @param {vec3} v0 
         * @param {vec3} v1 
         * @param {vec3} dst 
         * @returns vec3
         */
        subtract: function(v0, v1, dst){
            dst = dst || new Float32Array(3);
            
            dst[0] = v0[0] - v1[0];
            dst[1] = v0[1] - v1[1];
            dst[2] = v0[2] - v1[2];
            
            return dst;
        },
        /*
            i   j   k
            v00 v01 v02
            v10 v11 v12
        
        */
        /**
        * Calculates the cross product of two vec3s.
        * (v0 x v1)
        * 
        * @param {vec3} v0 
        * @param {vec3} v1 
        * @param {vec3} dst 
        * @returns vec3
        */
        cross: function(v0, v1, dst){
            dst = dst || new Float32Array(3);

            const v00 = v0[0];
            const v01 = v0[1];
            const v02 = v0[2];
            const v10 = v1[0];
            const v11 = v1[1];
            const v12 = v1[2];

            dst[0] = v01*v12 - v02*v11;
            dst[1] = -(v00*v12 - v02*v10);
            dst[2] = v00*v11 - v01*v10;
            
            return dst;
        }

    },
    
    vec4: {
        /**
         *  vec4 constructor
         * @param {array[4]} vals 
         * @returns vec4
         */
        construct: function(vals){
            vals = vals || [];
            const v4 = new Float32Array(3);

            v4[0] = vals[0] || 0;
            v4[1] = vals[1] || 0;
            v4[2] = vals[2] || 0;
            v4[3] = vals[3] || 0;

            return v4;
        },
    },

    /**
     *  Returns mat4 of all zeroes
     * @param {mat4} matrix 
     * @returns 
     */
    construct(matrix) {
        matrix = matrix || [];

        // init array to all 0s
        const m = new Float32Array(16);

        m[0] = matrix[0] || 0;
        m[1] = matrix[1] || 0;
        m[2] = matrix[2] || 0;
        m[3] = matrix[3] || 0;
        m[4] = matrix[4] || 0;
        m[5] = matrix[5] || 0;
        m[6] = matrix[6] || 0;
        m[7] = matrix[7] || 0;
        m[8] = matrix[8] || 0;
        m[9] = matrix[9] || 0;
        m[10] = matrix[10] || 0;
        m[11] = matrix[11] || 0;
        m[12] = matrix[12] || 0;
        m[13] = matrix[13] || 0;
        m[14] = matrix[14] || 0;
        m[15] = matrix[15] || 0;    

        return m;
    },

    /**
     * creates an mat4 indentity matrix
     * 
     * @returns mat4
     */
    identity: function(){
        const matrix = this.construct();
        
        // make pivots 1
        matrix[0] = 1;
        matrix[5] = 1;
        matrix[10] = 1;
        matrix[15] = 1;
        
        return matrix;
    },
    
    // calculates and returns inverse of given matrix
    //  does so by calculating the 1/determinant and adjoint
    /**
     * 
     * @param {mat4} matrix 
     * @param {mat4} invMat 
     * @returns 
     */
    inverse: function(matrix, invMat){
        // get vals
        const mx1 = matrix[0];
        const mx2 = matrix[1];
        const mx3 = matrix[2];
        const mx4 = matrix[3];
        
        const my1 = matrix[4];
        const my2 = matrix[5];
        const my3 = matrix[6];
        const my4 = matrix[7];
        
        const mz1 = matrix[8];
        const mz2 = matrix[9];
        const mz3 = matrix[10];
        const mz4 = matrix[11];
        
        const mw1 = matrix[12];
        const mw2 = matrix[13];
        const mw3 = matrix[14];
        const mw4 = matrix[15];
        
        // calc det A
        const tmp1 = (my2*(mz3*mw4-mw3*mz4)-mz2*(my3*mw4-mw3*my4)+mw2*(my3*mz4-mz3*my4));
        const tmp2 = (mx2*(mz3*mw4-mw3*mz4)-mz2*(mx3*mw4-mw3*mx4)+mw2*(mx3*mz4-mz3*mx4));
        const tmp3 = (mx2*(my3*mw4-mw3*my4)-my2*(mx3*mw4-mw3*mx4)+mw2*(mx3*my4-my3*mx4));
        const tmp4 = (mx2*(my3*mz4-mz3*my4)-my2*(mx3*mz4-mz3*mx4)+mz2*(mx3*my4-my3*mx4));
        
        const det = mx1*tmp1-my1*tmp2+mz1*tmp3-mw1*tmp4;
        
        // cofactor col 2
        const tmp5 = (my1*(mz3*mw4-mw3*mz4)-mz1*(my3*mw4-mw3*my4)+mw1*(my3*mz4-mz3*my4));
        const tmp6 = (mx1*(mz3*mw4-mw3*mz4)-mz1*(mx3*mw4-mw3*mx4)+mw1*(mx3*mz4-mz3*mx4));
        const tmp7 = (mx1*(my3*mw4-mw3*my4)-my1*(mx3*mw4-mw3*mx4)+mw1*(mx3*my4-my3*mx4));
        const tmp8 = (mx1*(my3*mz4-mz3*my4)-my1*(mx3*mz4-mz3*mx4)+mz1*(mx3*my4-my3*mx4));
        
        // cofactor col 3 (reverse sign since positions flipped)
        const tmp9  = -(my2*(mz1*mw4-mw1*mz4)-mz2*(my1*mw4-mw1*my4)+mw2*(my1*mz4-mz1*my4));
        const tmp10 = -(mx2*(mz1*mw4-mw1*mz4)-mz2*(mx1*mw4-mw1*mx4)+mw2*(mx1*mz4-mz1*mx4));
        const tmp11 = -(mx2*(my1*mw4-mw1*my4)-my2*(mx1*mw4-mw1*mx4)+mw2*(mx1*my4-my1*mx4));
        const tmp12 = -(mx2*(my1*mz4-mz1*my4)-my2*(mx1*mz4-mz1*mx4)+mz2*(mx1*my4-my1*mx4));
        
        // cofactor col 4
        const tmp13 = (my2*(mz3*mw1-mw3*mz1)-mz2*(my3*mw1-mw3*my1)+mw2*(my3*mz1-mz3*my1));
        const tmp14 = (mx2*(mz3*mw1-mw3*mz1)-mz2*(mx3*mw1-mw3*mx1)+mw2*(mx3*mz1-mz3*mx1));
        const tmp15 = (mx2*(my3*mw1-mw3*my1)-my2*(mx3*mw1-mw3*mx1)+mw2*(mx3*my1-my3*mx1));
        const tmp16 = (mx2*(my3*mz1-mz3*my1)-my2*(mx3*mz1-mz3*mx1)+mz2*(mx3*my1-my3*mx1));

        if(det == 0){return false}; // in case matrix is not invertible
        
        // get 1/d
        const dInv = 1.0 / det;
        invMat =  invMat || this.construct();
        
        // col 1
        invMat[0]  = tmp1 * dInv;
        invMat[1]  = -1 * tmp2 * dInv;
        invMat[2]  = tmp3 * dInv;
        invMat[3] = -1 * tmp4 * dInv;
        // col 2
        invMat[4]  = -1 * tmp5 * dInv;
        invMat[5]  = tmp6 * dInv;
        invMat[6]  = -1 * tmp7 * dInv;
        invMat[7] = tmp8 * dInv;
        // col 3
        invMat[8]  = tmp9 * dInv;
        invMat[9]  = -1 * tmp10 * dInv;
        invMat[10] = tmp11 * dInv;
        invMat[11] = -1 * tmp12 * dInv;
        // col 4
        invMat[12]  = -1 * tmp13 * dInv;
        invMat[13]  = tmp14 * dInv;
        invMat[14] = -1 * tmp15 * dInv;
        invMat[15] = tmp16 * dInv;
        
        return invMat;
    },
    
    /*
        Rotation about the x axis
    
    1 0    0     0
    0 cos@ -sin@ 0
    0 sin@ cos@  0
    0 0    0     1
    
    */
    /**
    * 
    * @param {mat4} matrix 
    * @param {float(radians)} rotation 
    * @param {mat4} dm 
    * @returns 
    */
    rotateX : function(matrix, rotation, dm){
        if(!rotation){return matrix};
        const c = Math.cos(rotation);
        const s = Math.sin(rotation);
        
        dm = dm || this.construct();
        
        dm[0] = matrix[0];
        dm[1] = matrix[1];
        dm[2] = matrix[2];
        dm[3] = matrix[3];
        
        dm[4] = matrix[4] * c + matrix[8] * s;
        dm[5] = matrix[5] * c + matrix[9] * s;
        dm[6] = matrix[6] * c + matrix[10] * s;
        dm[7] = matrix[7] * c + matrix[11] * s;
        
        dm[8] = matrix[4] * -s + matrix[8] * c;
        dm[9] = matrix[5] * -s + matrix[9] * c;
        dm[10] = matrix[6] * -s + matrix[10] * c;
        dm[11] = matrix[7] * -s + matrix[11] * c;
        
        dm[12] = matrix[12];
        dm[13] = matrix[13];
        dm[14] = matrix[14];
        dm[15] = matrix[15];
        
        return dm;
        
    },

    /*
        Rotation about the y axis

        cos@  0   sin@
        0     1   0
        -sin@ 0   cos@

    */
    /**
    * 
    * @param {mat4} matrix 
    * @param {float(radians)} rotation 
    * @param {mat4} dm 
    * @returns 
    */
    rotateY: function(matrix, rotation, dm){
        if(!rotation){return matrix};
        const c = Math.cos(rotation);
        const s = Math.sin(rotation);
        
        dm = dm || this.construct();
        
        dm[0] = matrix[0] * c - matrix[8] * s;
        dm[1] = matrix[1] * c - matrix[9] * s;
        dm[2] = matrix[2] * c - matrix[10] * s;
        dm[3] = matrix[3] * c - matrix[11] * s;
        
        dm[4] = matrix[4];
        dm[5] = matrix[5];
        dm[6] = matrix[6];
        dm[7] = matrix[7];
        
        dm[8] = matrix[0] * s + matrix[8] * c  ;
        dm[9] = matrix[1] * s + matrix[9] * c  ;
        dm[10] = matrix[2] * s + matrix[10] * c;
        dm[11] = matrix[3] * s + matrix[11] * c;
        
        dm[12] = matrix[12];
        dm[13] = matrix[13];
        dm[14] = matrix[14];
        dm[15] = matrix[15];
        
        return dm;
        
    },

    /*
        rotation about the z-axis

        cos@ -sin@ 0
        sin@ cos@  0
        0    0     1

    */
   /**
    * 
    * @param {mat4} matrix 
    * @param {float(radians)} rotation 
    * @param {mat4} dm 
    * @returns 
    */
    rotateZ: function(matrix, rotation, dm){
        if(!rotation){return matrix};
        const c = Math.cos(rotation);
        const s = Math.sin(rotation);
        
        dm = dm || this.construct();
        
        dm[0] = matrix[0] * c + matrix[4] * s;
        dm[1] = matrix[1] * c + matrix[5] * s;
        dm[2] = matrix[2] * c + matrix[6] * s;
        dm[3] = matrix[3] * c + matrix[7] * s;
        
        dm[4] = matrix[4] * c - matrix[0] * s;
        dm[5] = matrix[5] * c - matrix[1] * s;
        dm[6] = matrix[6] * c - matrix[2] * s;
        dm[7] = matrix[7] * c - matrix[3] * s;
        
        dm[8] =  matrix[8];
        dm[9] =  matrix[9];
        dm[10] = matrix[10];
        dm[11] = matrix[11];
        
        dm[12] = matrix[12];
        dm[13] = matrix[13];
        dm[14] = matrix[14];
        dm[15] = matrix[15];
        
        return dm;
    },

    /**
     * Calculates look at matrix for a set camera position and target
     * (inverse of view matrix)
     * 
     * @param {vec3} camPos 
     * @param {vec3} target 
     * @param {vec3} dst 
     * @returns 
     */
    lookAt: function(camPos, target, up, dst){
        dst = dst || new Float32Array(16);
        up = up || this.vec3.construct([0, 1, 0]);

        const fvec = this.vec3.normalize(this.vec3.subtract(camPos, target)); // forward vector
        const rvec = this.vec3.normalize(this.vec3.cross(fvec, up)); // right vector
        const uvec = this.vec3.normalize(this.vec3.cross(rvec, fvec)); // up vector

        dst[0] = rvec[0];
        dst[1] = rvec[1];
        dst[2] = rvec[2];
        dst[3] = 0;
        dst[4] = uvec[0];
        dst[5] = uvec[1];
        dst[6] = uvec[2];
        dst[7] = 0;
        dst[8] = fvec[0];
        dst[9] = fvec[1];
        dst[10] = fvec[2];
        dst[11] = 0;
        dst[12] = 0;
        dst[13] = 0;
        dst[14] = 0;
        dst[14] = 1;

        return dst;
    }
        
};

