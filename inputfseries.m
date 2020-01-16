% Find fourier series approximation of input signal

syms x f w %a b c h
a = .2;
b = 1.4*pi;
c = 2*pi;
h = 1;

func_period = 2*pi;
w0 = 2*pi/func_period;
plot_domain = 0:.01:2*pi;


f = h*triangularPulse(a,b,c,x);
f_FT = fourier(f, x, w);

%f_FT = subs(f_FT, a, 0);
%f_FT = subs(f_FT, b, 1);
%f_FT = subs(f_FT, c, .7);
%f_FT = subs(f_FT, h, 1);

plot_freqs = -50.05:.1:50.05;
f_ft_plot = subs(f_FT, w,plot_freqs);

figure(); 
subplot(2,1,1)
plot(plot_freqs, abs(f_ft_plot))
subplot(2,1,2)
plot(plot_freqs, angle(f_ft_plot))

num_cosine_terms = 20;

terms = [-num_cosine_terms:num_cosine_terms];
cn = zeros(length(terms),1);
for n = terms
    if n ~=0
        cn(n+num_cosine_terms+1) = (1/func_period) *...
                                   double(subs(f_FT,w,n*w0));
    else
        cn(n+num_cosine_terms+1) = (1/func_period) *...
                                   real(subs(f_FT,w,.00001));
    end
end

f_terms = cn(num_cosine_terms+1) * ones(1, length(plot_domain));
for index = 1:num_cosine_terms
  f_terms(index+1,:) = ...
    cn(end-index+1)*exp(1i*(num_cosine_terms - index+1)*plot_domain) + ...
    cn(index)*exp(-1i*(num_cosine_terms - index+1)*plot_domain);
end

% Convert to fourier series with fundamental frequency 2*pi
% w0 = 2*pi;
% cn = [subs(f_FT, w, .001)];
% cn = [subs(f_FT, w, [-num_cosine_terms:-1] * w0), cn,...
%       subs(f_FT, w, [1:num_cosine_terms] * w0)];
% cn = double(cn);

%f_terms = real(cn(num_cosine_terms+1)) * ones(1,length(plot_domain));
%for index = 1:num_cosine_terms
%    f_terms(index+1, :) = (cn(end-index)*exp(-1i*index*plot_domain) +...
%                              cn(index)*exp(1i*index*plot_domain));
%end
%f_terms = f_terms / w0;

figure(); 
subplot(2,1,1);
plot(sum(real(f_terms)))
hold on; plot(subs(f, plot_domain))
subplot(2,1,2);
plot(terms,real(cn))


