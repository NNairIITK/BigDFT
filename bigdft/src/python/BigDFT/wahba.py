from numpy import *
from math import sqrt

def apply_R(R,A):
    "Apply the rotation on the set of vectors A"
    A2 = R*A.T
    A2 = A2.T
    return A2

def apply_t(t,A):
    "Apply a translation on the set of vectors A"
    n=A.shape[0]
    return A + tile(t, (1,n)).T

def apply_Rt(R,t,A):
    "Rotate the element and apply the translation on the rotated vector"
    RA=apply_R(R,A)
    return apply_t(t,RA)

# Input: expects Nx3 matrix of points
# Returns R,t
# R = 3x3 rotation matrix
# t = 3x1 column vector

def rigid_transform_3D(A, B):
    "Find the transformation R and t such that R*A + t ~= B, with an error quantified by J"
    assert len(A) == len(B)

    N = A.shape[0]; # total points

    centroid_A = mean(A, axis=0)
    centroid_B = mean(B, axis=0)

    #print 'centre',centroid_A,centroid_B
    # centre the points
    AA = A - tile(centroid_A, (N, 1))
    BB = B - tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = transpose(AA) * BB

    #print 'H',H
    
    U, S, Vt = linalg.svd(H)

    R = Vt.T * U.T

    # special reflection case
    if linalg.det(R) < 0:
       print "#Reflection detected"
       Vt[2,:] *= -1
       R = Vt.T * U.T

    t = -R*centroid_A.T + centroid_B.T

    #print t

    #identify also the accuracy of wahba
    A2 = R*A.T + tile(t, (1, N))
    A2 = A2.T
    
    # Find the error
    err = A2 - B
    
    err = multiply(err, err)
    err = sum(err)
    rmse = sqrt(err/N);
    
    return R, t,rmse

# Test with random data
if __name__ == '__main__':
    # Random rotation and translation
    R = mat(random.rand(3,3))
    t = mat(random.rand(3,1))
    
    # make R a proper rotation matrix, force orthonormal
    U, S, Vt = linalg.svd(R)
    R = U*Vt
    
    # remove reflection
    if linalg.det(R) < 0:
        Vt[2,:] *= -1
        R = U*Vt
    
    # number of points
    n = 10
    
    A = mat(random.rand(n,3));
    B = R*A.T + tile(t, (1, n))
    B = B.T;

# recover the transformation
    ret_R, ret_t = rigid_transform_3D(A, B)

    A2 = (ret_R*A.T) + tile(ret_t, (1, n))
    A2 = A2.T
    
    # Find the error
    err = A2 - B
    
    err = multiply(err, err)
    err = sum(err)
    rmse = sqrt(err/n);
    
    print "Points A"
    print A
    print ""
    
    print "Points B"
    print B
    print ""
    
    print "Rotation"
    print R
    print ""
    
    print "Translation"
    print t
    print ""
    
    print "RMSE:", rmse
    print "If RMSE is near zero, the function is correct!"
    
