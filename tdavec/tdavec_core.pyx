#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
import numpy as np
cimport numpy as np
from cython.parallel import prange




def DiagToPD(D):
    """
    Generates a list of persistence diagrams (PD) from a given list of persistence diagrams (D).

    Parameters:
    - D (list): A list of persistence diagrams, where each diagram is represented as a numpy array.

    Returns:
    - PD (list): A list of persistence diagrams (PD), where each PD is represented as a numpy array.
      Each PD contains two columns: the first column represents the birth values of the persistence pairs,
      and the second column represents the persistence, i.e. the death values minus the birth values.
    """
    PD = [ np.transpose(np.array([D[dim][:,0], D[dim][:,1] - D[dim][:,0]])) for dim in range(len(D))]
    return PD


def computePersistenceBlock_dim0(np.ndarray[np.float64_t, ndim=1] x, 
                    np.ndarray[np.float64_t, ndim=1] y, 
                    np.ndarray[np.float64_t, ndim=1] ySeq, 
                    np.ndarray[np.float64_t, ndim=1] lam):
    """
    Compute the VPB values for dimension 0.

    Parameters:
        x (numpy.ndarray): The x values.
        y (numpy.ndarray): The y values.
        ySeq (numpy.ndarray): The sequence of y values.
        lam (numpy.ndarray): The lambda values.

    Returns:
        numpy.ndarray: The computed VPB values.
    """
    cdef Py_ssize_t n = ySeq.shape[0] - 1
    cdef int nPoints = y.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] vpb = np.zeros(n, dtype=np.float64)
    cdef int i, j
    cdef double c, d, y_cd, lam_cd, yMin, yMax
    for i in range(n):
        c = ySeq[i]
        d = ySeq[i+1]
        for j in range(nPoints):
            if c - lam[j] < y[j] < d + lam[j]:
                y_cd = y[j]
                lam_cd = lam[j]
                yMin = max(c, y_cd - lam_cd)
                yMax = min(d, y_cd + lam_cd)
                vpb[i] += 0.5 * (yMax**2 - yMin**2) / (ySeq[i+1] - ySeq[i])
    return vpb

def pmax(num, vec):
    """
    Compute the element-wise maximum of a scalar value and a NumPy array.

    Parameters:
        num (float): The scalar value.
        vec (numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: The resulting array with the element-wise maximum.
    """
    return np.array([max(num, vec[i_]) for i_ in range(vec.size)])

def pmin(num, vec):
    """
    Compute the element-wise minimum of a scalar value and a NumPy array.

    Parameters:
        num (float): The scalar value.
        vec (numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: The resulting array with the element-wise minimum.
    """
    return np.array([min(num, vec[i_]) for i_ in range(vec.size)])

def computePersistenceBlock_dim1(np.ndarray[np.float64_t, ndim=1] x, 
                    np.ndarray[np.float64_t, ndim=1] y, 
                    np.ndarray[np.float64_t, ndim=1] xSeq, 
                    np.ndarray[np.float64_t, ndim=1] ySeq, 
                    np.ndarray[np.float64_t, ndim=1] lam):
    """
    Compute the Vector Persistence Block (VPB) vectorization for a given set of points in dimension 1.

    Parameters:
        x (numpy.ndarray): The x-coordinates of the points.
        y (numpy.ndarray): The y-coordinates of the points.
        xSeq (numpy.ndarray): The x-coordinates of the grid points.
        ySeq (numpy.ndarray): The y-coordinates of the grid points.
        lam (numpy.ndarray): The lambda values.

    Returns:
        numpy.ndarray: The VPB matrix.
    """
    cdef int n = xSeq.shape[0] - 1
    cdef int m = ySeq.shape[0] - 1
    cdef np.ndarray[np.float64_t, ndim=1] dx = np.diff(xSeq)
    cdef np.ndarray[np.float64_t, ndim=1] dy = np.diff(ySeq)
    cdef np.ndarray[np.float64_t, ndim=2] vpb = np.zeros((n, m), dtype=np.float64)
    cdef int i, j, k, num_inds
    cdef double a, b, c, d, add
    cdef np.ndarray[np.int64_t, ndim=1] inds
    cdef np.ndarray[np.float64_t, ndim=1] xInd, yInd, lamInd, xMin, xMax, yMin, yMax

    for i in range(n):
        a, b = xSeq[i], xSeq[i+1]
        for j in range(m):
            c, d = ySeq[j], ySeq[j+1]
            # Using bitwise operations for logical conditions
            xCond = (x + lam >= a) & (x - lam <= b)
            yCond = (y + lam >= c) & (y - lam <= d)
            inds = np.where(xCond & yCond)[0]
            num_inds = inds.shape[0]

            if num_inds > 0:
                xInd = x.take(inds)
                yInd = y.take(inds)
                lamInd = lam.take(inds)
                xMin = np.maximum(a, xInd - lamInd)
                xMax = np.minimum(b, xInd + lamInd)
                yMin = np.maximum(c, yInd - lamInd)
                yMax = np.minimum(d, yInd + lamInd)
                add = 0.5 * np.sum((xMax - xMin) * (yMax - yMin) * (xMax + xMin + yMax + yMin)) / dx[i] / dy[j]
                vpb[i, j] += add
    return np.asarray(vpb)

def computePersistenceBlock(D, homDim, xSeq, ySeq, tau=0.3):
    """
    Compute the VPB vectorization using the given parameters.

    Parameters:
        D (list): Persistence Diagram (list of birth-death arrays for each dimension).
        homDim (int): The dimension along which the homogeneity is computed.
        xSeq (numpy.ndarray): The x-coordinates of the grid points.
        ySeq (numpy.ndarray): The y-coordinates of the grid points.
        tau (float, optional): The tau value. Defaults to 0.3.

    Returns:
        numpy.ndarray: The VPB matrix.
    """
    x = D[homDim][:,0]
    y = D[homDim][:,1] - x
    lam = tau * y
    if homDim == 0:
        return computePersistenceBlock_dim0(x, y, ySeq, lam)
    else:
        return computePersistenceBlock_dim1(x, y, xSeq, ySeq, lam)

def computePersistenceLandscape(D, homDim, scaleSeq, k=1):
    """
    Compute the persistence landscape (PL) for a given homological dimension, scale sequence, and order of landscape.

    Parameters:
        D (numpy.ndarray): Persistence Diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension along which the PL is computed.
        scaleSeq (numpy.ndarray): The sequence of scale values.
        k (int, optional): The order of the PL. Defaults to 1.

    Returns:
        numpy.ndarray: The persistence landscape vector.
    """
    birth, death = D[homDim][:,0], D[homDim][:,1]
    Lambda = [
        np.sort(pmax(0, np.apply_along_axis(min, 0, np.array([s-birth, death-s]))))[-k]
        for s in scaleSeq]
    return np.array(Lambda)

def computePersistenceSilhouette(D, homDim, scaleSeq, p=1):
    """
    Compute the Persistence Silhouette vectorization for a given homological dimension, scale sequence, and power.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension along which the PS is computed.
        scaleSeq (numpy.ndarray): The sequence of scale values.
        p (int, optional): The power to raise the difference between y and x. Defaults to 1.

    Returns:
        numpy.ndarray: The persistence spectrum vector.
    """
    x, y = D[homDim][:,0], D[homDim][:,1]
    pp = (y-x)**p
    w = pp/np.sum(pp)

    phi = []
    for k in range(len(scaleSeq)-1):
        alpha1 = pmax(scaleSeq[k], x)
        alpha2 = pmax(scaleSeq[k], (x+y)/2)
        beta1 = pmin(scaleSeq[k+1], (x+y)/2)
        beta2 = pmin(scaleSeq[k+1], y)
        b1 = pmax(0,beta1-alpha1)*((beta1+alpha1)/2-x)
        b2 = pmax(0,beta2-alpha2)*(y-(beta2+alpha2)/2)
        phi.append( np.sum(w*(b1+b2))/(scaleSeq[k+1]-scaleSeq[k]))
    return np.array(phi)

def computeNormalizedLife(D, homDim, scaleSeq):
    """
    Compute theNormalized Life Curve vectorization for a given homological dimension, scale sequence, and power.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension along which the NL is computed.
        scaleSeq (numpy.ndarray): The sequence of scale values.

    Returns:
        numpy.ndarray: The nonlinear landscape vector.
    """
    x, y = D[homDim][:,0], D[homDim][:,1]
    lL = (y-x)/sum(y-x)
    nl = []
    for k in range(len(scaleSeq)-1):
        b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
        nl.append( np.sum(lL*pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
    return np.array(nl)

def computeBettiCurve(D, homDim, scaleSeq):
    """
    Compute the Vector Summary of the Betti Curve    (VAB) vectorization for a given homological dimension, scale sequence, and power.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension along which the VAB is computed.
        scaleSeq (numpy.ndarray): The sequence of scale values.

    Returns:
        numpy.ndarray: The VAB vector.
    """
    x, y = D[homDim][:,0], D[homDim][:,1]
    vab = []
    for k in range( len(scaleSeq)-1):
        b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
        vab.append( sum(pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
    return np.array(vab)

def computeEulerCharacteristic(D, maxhomDim, scaleSeq):
    """
    Compute the Euler Characteristic Curve (ECC) vectorization for a given homological dimension, maximum homological dimension, and scale sequence.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        maxhomDim (int): The maximum homological dimension.
        scaleSeq (numpy.ndarray): The sequence of scale values.

    Returns:
        numpy.ndarray: The ECC vector.
    """
    ecc = np.zeros( len(scaleSeq)-1)
    for d in range(maxhomDim+1):
        ecc = ecc + (-1)**d * computeBettiCurve(D, d, scaleSeq)
    return ecc

def computePersistentEntropy(D, homDim, scaleSeq):
    """
    Compute the Persistence Entropy Summary (PES) vectorization for a given homological dimension, scale sequence, and persistence diagram.

    Parameters:
        D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
        homDim (int): The homological dimension.
        scaleSeq (numpy.ndarray): The sequence of scale values.

    Returns:
        list: The PES values.
    """
    x, y = D[homDim][:,0], D[homDim][:,1]
    lL = (y-x)/np.sum(y-x)
    entr = -lL*np.log10(lL)/np.log10(2)
    pes = []
    for k in range( len(scaleSeq)-1):
        b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
        pes.append( np.sum(entr*pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
    return pes

from scipy.stats import norm
def pnorm(x, mean, sd):
    """
    Calculate the cumulative distribution function of a normal distribution.

    Parameters:
        x (float): The value at which to calculate the cumulative distribution function.
        mean (float): The mean of the normal distribution.
        sd (float): The standard deviation of the normal distribution.

    Returns:
        float: The cumulative distribution function value at x.
    """
    return norm.cdf(x, mean, sd)

def outer(x, y):
    """
    Generate the outer product of two arrays.

    Parameters:
        x (array-like): The first input array.
        y (array-like): The second input array.

    Returns:
        numpy.ndarray: The outer product of the input arrays.
    """
    return np.array([x_*y_ for y_ in y for x_ in x])

def PSurfaceH0(point, y_lower, y_upper, sigma, maxP):
    """
    Calculate the surface probability density function for a specific homDim=0 point on the y-axis .

    Parameters:
        point (tuple): A tuple containing the x and y coordinates of the point.
        y_lower (float): The lower bound of the y-axis interval.
        y_upper (float): The upper bound of the y-axis interval.
        sigma (float): The standard deviation of the normal distribution.
        maxP (float): The maximum value of the y-axis.

    Returns:
        float: The surface probability density function value at the given point.

    """
    y = point[1]
    out2 = pnorm(y_upper, y, sigma) - pnorm(y_lower, y, sigma)
    wgt = y/maxP if y<maxP else 1
    return wgt*out2

def PSurfaceHk(point, y_lower, y_upper, x_lower, x_upper, sigma, maxP):
    """
    Calculate the surface probability density function for a specific homDim>0point in a two-dimensional space.

    Parameters:
        point (tuple): A tuple containing the x and y coordinates of the point.
        y_lower (float): The lower bound of the y-axis interval.
        y_upper (float): The upper bound of the y-axis interval.
        x_lower (float): The lower bound of the x-axis interval.
        x_upper (float): The upper bound of the x-axis interval.
        sigma (float): The standard deviation of the normal distribution.
        maxP (float): The maximum value of the y-axis.

    Returns:
        float: The surface probability density function value at the given point.

    """
    x, y = point[0], point[1]
    out1 = pnorm(x_upper,x,sigma) - pnorm(x_lower,x,sigma)
    out2 = pnorm(y_upper,y,sigma) - pnorm(y_lower,y,sigma)
    wgt = y/maxP if y<maxP else 1
    return wgt*outer(out1, out2)

def computePI(D, homDim, xSeq, ySeq, sigma):
    """
    Compute the surface Persistence Image (PI) vectorization for a given Persistence diagram

    Args:
        D (list): Persistence Diagram (list of birth-death arrays for each dimension).
        homDim (int): The dimension to compute the surface probability density function for.
        xSeq (list): The x-axis sequence.
        ySeq (list): The y-axis sequence.
        sigma (float): The standard deviation of the normal distribution.

    Returns:
        numpy.ndarray: The surface probability density function values for each data point.

    """
    D_ = D[homDim]
    D_[:,1] = D_[:,1] - D_[:,0]
    n_rows = D_.shape[0]

    resB = len(xSeq) - 1
    resP = len(ySeq)-1
    minP, maxP = ySeq[0], ySeq[-1]
    dy = (maxP-minP)/resP
    y_lower = np.arange(minP, maxP, dy)
    y_upper = y_lower + dy

    nSize = resP if homDim == 0 else resP*resB
    Psurf_mat = np.zeros( (nSize, n_rows))
    if homDim==0:
        for i in range(n_rows):
            Psurf_mat[:, i] = PSurfaceH0(D_[i, :], y_lower, y_upper, sigma, maxP)
    else:
        minB, maxB = xSeq[0], xSeq[-1]
        dx = (maxB-minB)/resB
        x_lower = np.arange(minB, maxB, dx)
        x_upper = x_lower + dx
        for i in range(n_rows):
            Psurf_mat[:, i] = PSurfaceHk(D_[i, :], y_lower, y_upper, x_lower, x_upper, sigma, maxP)
    out = np.sum(Psurf_mat, axis = 1)
    return out

def computeFDA(PD, maxD, homDim = 0, K = 10):
    X = np.zeros( (2*K+1))
    pd = PD[homDim]
    b = pd[:,0]/maxD; d = pd[:,1]/maxD
    X[0] = np.sum(d - b)
    for m in range(1, K+1):
        c = 2*m*np.pi
        alpha_sin = np.sin(c*d)-np.sin(c*b)
        alpha_cos = np.cos(c*d)-np.cos(c*b)
        X[2*m-1] = -np.sqrt(2)/c * np.sum(alpha_cos)
        X[2*m] = np.sqrt(2)/c * np.sum(alpha_sin)
    return X