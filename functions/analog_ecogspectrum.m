function [ecog_spectrum] = analog_ecogspectrum(samplingfreq,ECoG)
% ECoG processing
    detrended_ECoG = detrend(ECoG,'constant');
    ecog_spectrum.params.tapers = [5,9];
    ecog_spectrum.params.Fs = samplingfreq;
    ecog_spectrum.params.fpass = [1,100];
    ecog_spectrum.movingwin = [5,1/5];
    [S5,T5,F5] = mtspecgramc(detrended_ECoG,ecog_spectrum.movingwin,ecog_spectrum.params);
    % Output:
    %       S       (spectrum in form time x frequency x channels/trials if trialave=0; 
    %               in the form time x frequency if trialave=1)
    %       t       (times)
    %       f       (frequencies)
    %       Serr    (error bars) only for err(1)>=1
    
    normS5 = (S5-min(S5,[],'all'))./(max(S5,[],'all')-min(S5,[],'all'));
    % output
    ecog_spectrum.log_norm_spectrum = log(normS5)';
    ecog_spectrum.T = T5;
    ecog_spectrum.F = F5;
    ecog_spectrum.S = S5;
end

function [S,t,f,Serr] = mtspecgramc(data,movingwin,params)
    %________________________________________________________________________________________________________________________
    %%
    % Multi-taper time-frequency spectrum - continuous process
    %
    % Usage:
    % [S,t,f,Serr]=mtspecgramc(data,movingwin,params)
    % Input: 
    % Note units have to be consistent. Thus, if movingwin is in seconds, Fs
    % has to be in Hz. see chronux.m for more information.
    %       data        (in form samples x channels/trials) -- required
    %       movingwin         (in the form [window winstep] i.e length of moving
    %                                                 window and step size)
    %                                                 Note that units here have
    %                                                 to be consistent with
    %                                                 units of Fs - required
    %       params: structure with fields tapers, pad, Fs, fpass, err, trialave
    %       - optional
    %           tapers : precalculated tapers from dpss or in the one of the following
    %                    forms: 
    %                    (1) A numeric vector [TW K] where TW is the
    %                        time-bandwidth product and K is the number of
    %                        tapers to be used (less than or equal to
    %                        2TW-1). 
    %                    (2) A numeric vector [W T p] where W is the
    %                        bandwidth, T is the duration of the data and p 
    %                        is an integer such that 2TW-p tapers are used. In
    %                        this form there is no default i.e. to specify
    %                        the bandwidth, you have to specify T and p as
    %                        well. Note that the units of W and T have to be
    %                        consistent: if W is in Hz, T must be in seconds
    %                        and vice versa. Note that these units must also
    %                        be consistent with the units of params.Fs: W can
    %                        be in Hz if and only if params.Fs is in Hz.
    %                        The default is to use form 1 with TW=3 and K=5
    %                     Note that T has to be equal to movingwin(1).
    %
    %	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
    %                    -1 corresponds to no padding, 0 corresponds to padding
    %                    to the next highest power of 2 etc.
    %			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
    %			      	 to 512 points, if pad=1, we pad to 1024 points etc.
    %			      	 Defaults to 0.
    %           Fs   (sampling frequency) - optional. Default 1.
    %           fpass    (frequency band to be used in the calculation in the form
    %                                   [fmin fmax])- optional. 
    %                                   Default all frequencies between 0 and Fs/2
    %           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
    %                                   [0 p] or 0 - no error bars) - optional. Default 0.
    %           trialave (average over trials/channels when 1, don't average when 0) - optional. Default 0
    % Output:
    %       S       (spectrum in form time x frequency x channels/trials if trialave=0; 
    %               in the form time x frequency if trialave=1)
    %       t       (times)
    %       f       (frequencies)
    %       Serr    (error bars) only for err(1)>=1
    
    if nargin < 2; error('Need data and window parameters'); end;
    if nargin < 3; params=[]; end;
    
    if length(params.tapers)==3 & movingwin(1)~=params.tapers(2);
        error('Duration of data in params.tapers is inconsistent with movingwin(1), modify params.tapers(2) to proceed')
    end
    
    [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
    if nargout > 3 && err(1)==0; 
    %   Cannot compute error bars with err(1)=0. change params and run again.
        error('When Serr is desired, err(1) has to be non-zero.');
    end;
    data=change_row_to_column(data);
    [N,Ch]=size(data);
    Nwin=round(Fs*movingwin(1)); % number of samples in window
    Nstep=round(movingwin(2)*Fs); % number of samples to step through
    nfft=max(2.^(nextpow2(Nwin)+pad),Nwin);
    f=getfgrid(Fs,nfft,fpass); Nf=length(f);
    params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers
    
    winstart=1:Nstep:N-Nwin+1;
    nw=length(winstart); 
    
    if trialave
        S = zeros(nw,Nf);
        if nargout==4; Serr=zeros(2,nw,Nf); end;
    else
        S = zeros(nw,Nf,Ch);
        if nargout==4; Serr=zeros(2,nw,Nf,Ch); end;
    end
    
    for n=1:nw;
       indx=winstart(n):winstart(n)+Nwin-1;
       datawin=data(indx,:);
       if nargout==4
         [s,f,serr]=mtspectrumc(datawin,params);
         Serr(1,n,:,:)=squeeze(serr(1,:,:));
         Serr(2,n,:,:)=squeeze(serr(2,:,:));
       else
         [s,f]=mtspectrumc(datawin,params);
       end
       S(n,:,:)=s;
    end;
    S=squeeze(S); 
    if nargout==4;Serr=squeeze(Serr);end;
    winmid=winstart+round(Nwin/2);
    t=winmid/Fs;
end
function [S,f,Serr] = mtspectrumc(data,params)
    
    % Multi-taper spectrum - continuous process
    %
    % Usage:
    %
    % [S,f,Serr]=mtspectrumc(data,params)
    % Input: 
    % Note units have to be consistent. See chronux.m for more information.
    %       data (in form samples x channels/trials) -- required
    %       params: structure with fields tapers, pad, Fs, fpass, err, trialave
    %       -optional
    %           tapers : precalculated tapers from dpss or in the one of the following
    %                    forms: 
    %                    (1) A numeric vector [TW K] where TW is the
    %                        time-bandwidth product and K is the number of
    %                        tapers to be used (less than or equal to
    %                        2TW-1). 
    %                    (2) A numeric vector [W T p] where W is the
    %                        bandwidth, T is the duration of the data and p 
    %                        is an integer such that 2TW-p tapers are used. In
    %                        this form there is no default i.e. to specify
    %                        the bandwidth, you have to specify T and p as
    %                        well. Note that the units of W and T have to be
    %                        consistent: if W is in Hz, T must be in seconds
    %                        and vice versa. Note that these units must also
    %                        be consistent with the units of params.Fs: W can
    %                        be in Hz if and only if params.Fs is in Hz.
    %                        The default is to use form 1 with TW=3 and K=5
    %
    %	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
    %                    -1 corresponds to no padding, 0 corresponds to padding
    %                    to the next highest power of 2 etc.
    %			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
    %			      	 to 512 points, if pad=1, we pad to 1024 points etc.
    %			      	 Defaults to 0.
    %           Fs   (sampling frequency) - optional. Default 1.
    %           fpass    (frequency band to be used in the calculation in the form
    %                                   [fmin fmax])- optional. 
    %                                   Default all frequencies between 0 and Fs/2
    %           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
    %                                   [0 p] or 0 - no error bars) - optional. Default 0.
    %           trialave (average over trials/channels when 1, don't average when 0) - optional. Default 0
    % Output:
    %       S       (spectrum in form frequency x channels/trials if trialave=0; 
    %               in the form frequency if trialave=1)
    %       f       (frequencies)
    %       Serr    (error bars) only for err(1)>=1
    
    if nargin < 1; error('Need data'); end;
    if nargin < 2; params=[]; end;
    [tapers,pad,Fs,fpass,err,trialave,params] = getparams(params);
    if nargout > 2 && err(1)==0; 
    %   Cannot compute error bars with err(1)=0. Change params and run again. 
        error('When Serr is desired, err(1) has to be non-zero.');
    end;
    data=change_row_to_column(data);
    N=size(data,1);
    nfft=max(2^(nextpow2(N)+pad),N);
    [f,findx]=getfgrid(Fs,nfft,fpass); 
    tapers=dpsschk(tapers,N,Fs); % check tapers
    J=mtfftc(data,tapers,nfft,Fs);
    J=J(findx,:,:);
    S=squeeze(mean(conj(J).*J,2));
    if trialave; S=squeeze(mean(S,2));end;
    if nargout==3; 
       Serr=specerr(S,J,err,trialave);
    end;
end
function J = mtfftc(data,tapers,nfft,Fs)

% Multi-taper fourier transform - continuous data
%
% Usage:
% J=mtfftc(data,tapers,nfft,Fs) - all arguments required
% Input: 
%       data (in form samples x channels/trials or a single vector) 
%       tapers (precalculated tapers from dpss) 
%       nfft (length of padded data)
%       Fs   (sampling frequency)
%                                   
% Output:
%       J (fft in form frequency index x taper index x channels/trials)
if nargin < 4; error('Need all input arguments'); end;
data=change_row_to_column(data);
[NC,C]=size(data); % size of data
[NK K]=size(tapers); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
data=data(:,:,ones(1,K)); % add taper indices to data
data=permute(data,[1 3 2]); % reshape data to get dimensions to match those of tapers
data_proj=data.*tapers; % product of data with tapers
J=fft(data_proj,nfft)/Fs;   % fft of projected data
end

function [tapers,pad,Fs,fpass,err,trialave,params] = getparams(params)

% Helper function to convert structure params to variables used by the
% various routines - also performs checks to ensure that parameters are
% defined; returns default values if they are not defined.
%
% Usage: [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params)
%
% Inputs:
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%           - optional
%             tapers : precalculated tapers from dpss or in the one of the following
%                       forms:  
%                       (1) A numeric vector [TW K] where TW is the
%                           time-bandwidth product and K is the number of
%                           tapers to be used (less than or equal to
%                           2TW-1). 
%                       (2) A numeric vector [W T p] where W is the
%                           bandwidth, T is the duration of the data and p 
%                           is an integer such that 2TW-p tapers are used. In
%                           this form there is no default i.e. to specify
%                           the bandwidth, you have to specify T and p as
%                           well. Note that the units of W and T have to be
%			                consistent: if W is in Hz, T must be in seconds
% 			                and vice versa. Note that these units must also
%			                be consistent with the units of params.Fs: W can
%		    	            be in Hz if and only if params.Fs is in Hz.
%                           The default is to use form 1 with TW=3 and K=5
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Outputs: 
% The fields listed above as well as the struct params. The fields are used
% by some routines and the struct is used by others. Though returning both
% involves overhead, it is a safer, simpler thing to do.

if ~isfield(params,'tapers') || isempty(params.tapers);  %If the tapers don't exist
     display('tapers unspecified, defaulting to params.tapers=[3 5]');
     params.tapers=[3 5];
end;
if ~isempty(params) && length(params.tapers)==3 
    % Compute timebandwidth product
    TW = params.tapers(2)*params.tapers(1);
    % Compute number of tapers
    K  = floor(2*TW - params.tapers(3));
    params.tapers = [TW  K];
end

if ~isfield(params,'pad') || isempty(params.pad);
    params.pad=0;
end;
if ~isfield(params,'Fs') || isempty(params.Fs);
    params.Fs=1;
end;
if ~isfield(params,'fpass') || isempty(params.fpass);
    params.fpass=[0 params.Fs/2];
end;
if ~isfield(params,'err') || isempty(params.err);
    params.err=0;
end;
if ~isfield(params,'trialave') || isempty(params.trialave);
    params.trialave=0;
end;

tapers=params.tapers;
pad=params.pad;
Fs=params.Fs;
fpass=params.fpass;
err=params.err;
trialave=params.trialave;
end
function data = change_row_to_column(data)

    % Helper routine to transform 1d arrays into column vectors that are needed
    % by other routines in Chronux
    %
    % Usage: data=change_row_to_column(data)
    % 
    % Inputs:
    % data -- required. If data is a matrix, it is assumed that it is of the
    % form samples x channels/trials and it is returned without change. If it
    % is a vector, it is transformed to a column vector. If it is a struct
    % array of dimension 1, it is again returned as a column vector. If it is a
    % struct array with multiple dimensions, it is returned without change
    % Note that the routine only looks at the first field of a struct array.
    % 
    % Ouputs:
    % data (in the form samples x channels/trials)
    %
    dtmp=[];
    if isstruct(data);
       C=length(data);
       if C==1;
          fnames=fieldnames(data);
          eval(['dtmp=data.' fnames{1} ';'])
          data=dtmp(:);
       end
    else
      [N,C]=size(data);
      if N==1 || C==1;
        data=data(:);
      end;
    end;
end

function [f,findx] = getfgrid(Fs,nfft,fpass)

% Helper function that gets the frequency grid associated with a given fft based computation
% Called by spectral estimation routines to generate the frequency axes 
% Usage: [f,findx]=getfgrid(Fs,nfft,fpass)
% Inputs:
% Fs        (sampling frequency associated with the data)-required
% nfft      (number of points in fft)-required
% fpass     (band of frequencies at which the fft is being calculated [fmin fmax] in Hz)-required
% Outputs:
% f         (frequencies)
% findx     (index of the frequencies in the full frequency grid). e.g.: If
% Fs=1000, and nfft=1048, an fft calculation generates 512 frequencies
% between 0 and 500 (i.e. Fs/2) Hz. Now if fpass=[0 100], findx will
% contain the indices in the frequency grid corresponding to frequencies <
% 100 Hz. In the case fpass=[0 500], findx=[1 512].
if nargin < 3; error('Need all arguments'); end;
df=Fs/nfft;
f=0:df:Fs; % all possible frequencies
f=f(1:nfft);
if length(fpass)~=1;
   findx=find(f>=fpass(1) & f<=fpass(end));
else
   [fmin,findx]=min(abs(f-fpass));
   clear fmin
end;
f=f(findx);

end

function [tapers,eigs] = dpsschk(tapers,N,Fs)

% Helper function to calculate tapers and, if precalculated tapers are supplied, 
% to check that they (the precalculated tapers) the same length in time as
% the time series being studied. The length of the time series is specified
% as the second input argument N. Thus if precalculated tapers have
% dimensions [N1 K], we require that N1=N.
% Usage: tapers=dpsschk(tapers,N,Fs)
% Inputs:
% tapers        (tapers in the form of: 
%                                   (i) precalculated tapers or,
%                                   (ii) [NW K] - time-bandwidth product, number of tapers) 
%
% N             (number of samples)
% Fs            (sampling frequency - this is required for nomalization of
%                                     tapers: we need tapers to be such
%                                     that integral of the square of each taper equals 1
%                                     dpss computes tapers such that the
%                                     SUM of squares equals 1 - so we need
%                                     to multiply the dpss computed tapers
%                                     by sqrt(Fs) to get the right
%                                     normalization)
% Outputs: 
% tapers        (calculated or precalculated tapers)
% eigs          (eigenvalues) 
if nargin < 3; error('Need all arguments'); end
sz=size(tapers);
if sz(1)==1 && sz(2)==2;
    [tapers,eigs]=dpss(N,tapers(1),tapers(2));
    tapers = tapers*sqrt(Fs);
elseif N~=sz(1);
    error('seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
end;
end

