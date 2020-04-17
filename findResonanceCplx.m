%% Find the complex resonance frequencies using Mike's method
function omega0 = findResonanceCplx(L,H,th,guess)

envvar = [pwd,'/include'];
setenv('FF_INCLUDEPATH',envvar);
file = 'iceSplineNoSolve.edp';
ffpp=['/usr/local/ff++/openmpi-2.1/3.61-1/bin/FreeFem++ -nw -ne ', file];

omega0 = guess; %Initial Guesses.

dw = 1e-2;
tol = 1e-6;
while (abs(dw)>tol && ~isinf(abs(dw)))
    %% First omega;
    T1 = 2*pi/(omega0);
    cmd=[ffpp,' -Tr ',num2str(real(T1)),' -Ti ',num2str(imag(T1)),...
        ' -H ',num2str(H), ' -L ',num2str(L),' -h ',num2str(th),' -N ',...
        num2str(3),' -isUni ',num2str(1)];
    [aa,bb1]=system(cmd);
    if(aa)
        disp(bb1)
        error('Cannot run program. Check path of FF++ or install it');
    end
    HReal = load('1_NoSolve/2_ModesMatrix/ReH0.dat');
    HImag = load('1_NoSolve/2_ModesMatrix/ImH0.dat');
    Hw = HReal + 1i*HImag;       
    
    %% Second omega
    T2 = 2*pi/(omega0+dw);
    cmd=[ffpp,' -Tr ',num2str(real(T2)),' -Ti ',num2str(imag(T2)),...
        ' -H ',num2str(H), ' -L ',num2str(L),' -h ',num2str(th),' -N ',...
        num2str(3),' -isUni ',num2str(1)];
    [aa,bb2]=system(cmd);
    if(aa)
        disp(bb2)
        error('Cannot run program. Check path of FF++ or install it');
    end
    HReal = load('1_NoSolve/2_ModesMatrix/ReH0.dat');
    HImag = load('1_NoSolve/2_ModesMatrix/ImH0.dat');
    Hw1 = HReal + 1i*HImag;
    
    %% Find the update.
    dH = (Hw1-Hw)/(dw);
    domega = eig(Hw,-dH);    
    dw = min(domega);
    if(~isinf(abs(dw)))
        omega0 = omega0 + dw;
    end
end

end