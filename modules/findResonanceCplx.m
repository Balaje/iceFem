%% Find the complex resonance frequencies using Mike's method
function omega0 = findResonanceCplx(L,H,th,guess,file,SolDir)
global ff

envvar = [pwd,'/include'];
setenv('FF_INCLUDEPATH',envvar);
ffpp=[ff,' -nw -ne iceSplineNoSolve.edp'];

omega0 = guess; %Initial Guesses.

dw = 1e-2;
tol = 1e-6;
while (abs(dw)>tol && ~isinf(abs(dw)))
    %% First omega;
    T1 = 2*pi/(omega0);
    cmd=[ffpp,' -Tr ',num2str(real(T1)),' -Ti ',num2str(imag(T1)),...
        ' -H ',num2str(H), ' -L ',num2str(L),' -h ',num2str(th),' -N ',...
        num2str(3),' -isUniIce ',num2str(1),' -isUniCav ',num2str(1)];
    [aa,bb1]=system(cmd);
    if(aa)
        disp(bb1)
        error('Cannot run program. Check path of FF++ or install it');
    end
    HReal = load([SolDir,'ReH0.dat']);
    HImag = load([SolDir,'ImH0.dat']);
    Hw = HReal + 1i*HImag;       
    
    %% Second omega
    T2 = 2*pi/(omega0+dw);
    cmd=[ffpp,' -Tr ',num2str(real(T2)),' -Ti ',num2str(imag(T2)),...
        ' -H ',num2str(H), ' -L ',num2str(L),' -h ',num2str(th),' -N ',...
        num2str(3),' -isUniIce ',num2str(1),' -isUniCav ',num2str(1)];
    [aa,bb2]=system(cmd);
    if(aa)
        disp(bb2)
        error('Cannot run program. Check path of FF++ or install it');
    end    
    HReal = load([SolDir,'ReH0.dat']);
    HImag = load([SolDir,'ImH0.dat']);
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