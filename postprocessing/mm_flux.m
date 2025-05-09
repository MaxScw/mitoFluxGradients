function flux = mm_flux(data, param, c)
  vmax = param(1, :);
  km = param(2, :);
  flux = (vmax.*data)./(data + km);
end