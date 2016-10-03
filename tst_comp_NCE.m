
aNCE_GS = zeros(size(time.tspan));
aNCE_AGS = zeros(size(time.tspan));

for i = 1 : mcRuns
    
    NCE_GS = -cell2mat(saveALL.F2{i}.ll)./pf.no_particles;
    NCE_AGS = -cell2mat(saveALL.F4{i}.ll)./pf.no_particles;
    norm_tmp = sqrt(NCE_GS.^2 + NCE_AGS.^2);
    
    aNCE_GS = aNCE_GS + NCE_GS./norm_tmp;
    aNCE_AGS = aNCE_AGS + NCE_AGS./norm_tmp;
    
    plot(time.tspan,NCE_GS,'b',time.tspan,NCE_AGS,'r');
    pause
    
end;