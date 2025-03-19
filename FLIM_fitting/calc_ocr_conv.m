temps = [22 28 31 36 37 40];
ocr = [0.5475 0.7908 1.2110 2.3672 2.68 2.9028]; % fmol/s

function k_ubs=get_k_ub(ocr, V)

k_ubs = 10^5.*2.*(ocr)./V;

end

V = 9.5*10^4; %muM
k_ubs = get_k_ub(ocr,V);

for k_ind=1:numel(k_ubs)
    
    name_str = "ocr_conv_T"+string(temps(k_ind))+"C.mat";
    k_ub = k_ubs(k_ind);
    save(name_str, "k_ub")

end