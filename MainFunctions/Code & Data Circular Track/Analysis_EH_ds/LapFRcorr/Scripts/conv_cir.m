function c = conv_cir(a,b)

%% Inputs
% a is the circular variable. N X 1
% b is the kernel. N X 1. Should be odd number

if size(a,2)>size(a,1)
    a = a';
end
if size(b,2)>size(b,1)
    b = b';
end

window = floor(length(b)/2);
c = zeros(size(a));
for ii = 1:length(a)
    range = ii-window:ii+window;
    range(range<1) = length(a)+range(range<1);
    range(range>length(a)) = range(range>length(a))-length(a);
    
    c(ii) = a(range)'*b;
end
    