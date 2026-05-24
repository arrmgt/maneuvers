function y=change_rate(x,irate,orate)
x=double(x);
if(irate==orate | length(x)==1)
    y=x;
else
    [minterp,mdecim]=interp_decim(irate,orate);
    blurf=Ninterp(x,minterp);
    y=decimateByFactors(blurf,mdecim,'FIR');
end
return