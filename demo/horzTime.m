function [timevector,varargout] = horzTime(poolsize,Z,R)
% [timevector [,A,H,D]] = horzTime(poolsize,Z,R)
%time the horizon program for a geographic grid

rotang = 0;
timevector = table;
for k=1:1+length(poolsize)
    p = k-1;
    if k==1
        tic
        [SF,SB] = horizonRotateLatLon(rotang,Z,R,false);
        elapsedTime = toc;
        timevector = [timevector; table(1,elapsedTime,"oneAz_NoP")]; %#ok<*AGROW>
        nout = max(nargout,1) - 1;
        if nout>=1
            S.Az = A;
            S.Horz = H;
            S.Dis = D;
            S.SForward = SF;
            S.SBackward = SB;
            varargout{1} = S;
        end
    else
        parpool(poolsize(p))
        tic
        if contains(class(R),'Geographic','IgnoreCase',true)
            [SF,SB] = horizonRotateLatLon(rotang,Z,R,true); %#ok<ASGLU>
        else
            [SF,SB] = horizonRotateProj(rotang,Z,R,true); %#ok<ASGLU>
        end
        elapsedTime = toc;
        timevector = [timevector; table(poolsize(p),elapsedTime,"oneAz")];
    end
    delete(gcp('nocreate'));
end
timevector.Properties.VariableNames = {'poolsize','elapsedTime','parallel'};

end