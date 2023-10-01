function NL = noiseDensity(f,s,w)

Nt = 17 - 30*log10(f); % terbulance 
Ns = 40 + 20*(s - 0.5) + 26*log10(f) - 60*log10( f + 0.03); % shipping
Nw = 50 + 7.5*sqrt(w) + 20*log10(f) - 40*log10(f + 0.4); % wind
Nth = -15 + 20*log10(f); % thermal

NL = 10*log10(10.^(Nt/10) + 10.^(Ns/10) + 10.^(Nw/10) + 10.^(Nth/10));
end

