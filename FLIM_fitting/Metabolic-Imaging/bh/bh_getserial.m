function y = getserial(filename)
    if ~bh_isvalid(filename)
        y = [];
        return
    end
    sp = bh_readsetup(filename);
    n = sp.no_of_meas_desc_blocks;
    m = sp.modules;
    
    sn = cell(1, n);
    for i=1:n
        m = bh_getmeasdesc(sp,i);
        sn{1,i} = m.mod_ser_no;
    end
    y = unique(sn);
    
end