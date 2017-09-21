function interpolated = pureGridder1d(zSamples,density,n)

% data must be in [-0.5, 0.5]
alpha=2; W=5; error=1e-3; S=ceil(0.91/error/alpha); beta=pi*sqrt((W/alpha*(alpha-1/2))^2-0.8);
% generate interpolating low-pass filter
Gz=alpha*n; F_kbZ=besseli(0,beta*sqrt(1-([-1:2/(S*W):1]).^2)); z=-(alpha*n/2):(alpha*n/2)-1; F_kbZ=(Gz/W)*besseli(0,beta*sqrt(1-([-1:2/(S*W):1]).^2));
%F_kbZ=F_kbZ/max(F_kbZ(:));
% generate Fourier transform of 1d interpolating kernels
kbZ=sqrt(alpha*n)*sin(sqrt((pi*W*z/Gz).^2-beta^2))./sqrt((pi*W*z/Gz).^2-beta^2);
% zero out output array in the alphaX grid - use a vector representation to be able to use sparse matrix structure
n = alpha*n; interpolated = zeros(n,1);
% convert samples to matrix indices
nz = (n/2+1) + n*zSamples;
% loop over samples in kernel
for lz = -(W-1)/2:(W-1)/2,
      % find nearest samples
      nzt = round(nz+lz);
      % compute weighting for KB kernel using linear interpolation
      zpos=S*((nz-nzt)-(-W/2))+1; kwz = F_kbZ(round(zpos))';
      % map samples outside the matrix to the edges
      nzt = max(nzt,1); nzt = min(nzt,n);
      % use sparse matrix to turn k-space trajectory into 2D matrix
      interpolated = interpolated+sparse(nzt,1,density.*kwz,n,1);
end;
interpolated = full(interpolated);
% zero out edge samples, since these may be due to samples outside the matrix
interpolated(1) = 0; interpolated(n) = 0;
n = n/alpha;
interpolated = myifft(interpolated,n);
% generate Fourier domain deapodization function and deapodize
deapod = kbZ(n/2+1:3*n/2)'; interpolated = interpolated./deapod;
% return to spatial domain
interpolated = abs(myfft3(interpolated));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=myifft(F,u)
% does the centering(shifts) and normalization
%needs a slight modification to actually support rectangular images
if nargin<2
f=ifftshift(ifft(ifftshift(F)))*sqrt(numel(F));
else
f=ifftshift(ifft(ifftshift(F)))*sqrt(u);
zoffset=ceil((size(f,1)-u)/2);
f=f(zoffset+1:zoffset+u);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=myfft3(f,u,v,w)
% does the centering(shifts) and normalization
%needs a slight modification to actually support rectangular images
if nargin<4
F=fftshift(fftn(fftshift(f)))/sqrt(prod(size(f)));
else
F=zeros(u,v,w);
zoffset=(u-size(f,1))/2; xoffset=(v-size(f,2))/2; yoffset=(w-size(f,3))/2;
F(zoffset+1:zoffset+size(f,1),xoffset+1:xoffset+size(f,2),yoffset+1:yoffset+size(f,3))=f;
F=fftshift(fftn(fftshift(F)))/sqrt(prod(size(f)));
end
