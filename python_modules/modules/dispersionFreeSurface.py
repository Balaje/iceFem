import numpy as np
import cmath

pi = np.pi

def oneroot(alpha, guess):
    ans1 = guess+1
    out = guess
    while(abs(ans1-out) > 1e-9):
        ans1 = out
        f = ans1*np.tanh(ans1) - alpha
        difff = np.tanh(ans1) + ans1*(1/np.cosh(ans1))**2
        out = ans1 - f/difff
    return out

def genalphastep(start, last, step):
    length = int(np.floor((last-start)/step))
    alphastep = np.zeros((length+2,1),dtype=complex)
    for m in np.arange(0,length+2):
        alphastep[m] = start + m*step
    alphastep[length+1] = abs(last)
    return alphastep

def homotopy(alpha, N1):
    if(N1 == 0):
        mroot = oneroot(1,1)
    else:
        mroot = oneroot(1, 1j*N1*pi)
    step = 0.043
    alphastep = genalphastep(1, abs(alpha), -step*(abs(alpha)<1) + step*(abs(alpha)>=1))
    for m in np.arange(1,len(alphastep)):
        mroot = oneroot(alphastep[m], mroot)
    argstep = genalphastep(0, cmath.phase(alpha), -(pi/30)*(cmath.phase(alpha)<0) + (pi/30)*(cmath.phase(alpha)>=0))
    argstep[len(argstep) - 1] = cmath.phase(alpha)
    newalphastep = np.zeros((len(argstep),1),dtype=complex)
    for m in np.arange(1, len(argstep)):
        newalphastep[m] = abs(alpha)*np.exp(1j*argstep[m])
        mroots = oneroot(newalphastep[m], mroot)
    return mroot

# The main dispersion equation function
def dispersion_free_surface(alpha, N, h):
    alpha = h*alpha
    mroots = np.zeros((N+1,1), dtype=complex)
    if(N==0):
        count = 0
        mroots[count] = homotopy(alpha, count)
    else:
        count = 0
        mroots[count] = homotopy(alpha, count)
        count = count + 1
        while(0<=1):
            mroots[count] = homotopy(alpha, count)
            if(abs(mroots[count] - (1j*count*pi + alpha/(1j*count*pi))) < 0.01):
                while(0<=1):
                    mroots[count] = oneroot(alpha, 1j*count*pi + alpha/(1j*count*pi))
                    if(abs(mroots[count] - (1j*count*pi + alpha/(1j*count*pi))) < 1e-15):
                        for m in np.arange(0,N-count):
                            mroots[count+m+1] = 1j*(count+m)*pi + alpha/(1j*(count+m)*pi)
                        count = N
                        break
                    if(count==N):
                        break
                    count = count+1
            if(count == N):
                break
            count = count + 1

        mroots = -(1j/h)*mroots
        # mroots[0] = -mroots[0]
        return mroots
