function alpha=waterAtten(f,temp,salinity,depth,pH)
    % Attenuation in fresh and sea water (From P.226 in Kinsler and Frey 4th)

    f_1=780*exp(temp/29); % Relaxation frequency for boric acid
    f_2=42e3*exp(temp/18); % Relaxation frequency for MgSO4
    A_a=83e-6*(salinity/35)*exp(temp/31-depth/91/1000+1.8*(pH-8)); % boric acid attenuation factor
    B_a=22e-3*(salinity/35)*exp(temp/14-depth/6/1000); % MgSO4 attenuation factor
    C_a=4.9e-13*exp(-temp/26-depth/25); % Pure water attenuation coefficient
    
    
    alpha=(A_a./(f_1^2+f.^2)+B_a./(f_2^2+f.^2)+C_a).*f.^2; %dB/m
end