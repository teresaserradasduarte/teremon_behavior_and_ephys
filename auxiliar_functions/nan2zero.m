function yy = nan2zero(xx)
yy=xx;
yy(isnan(yy))=0;
end

