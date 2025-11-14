function [par,nunique]=refl_voro_to_lay(voro,voroidx,k,NPV,NVMX);
  voroz  = zeros(NVMX,NPV,2);
  niface = zeros(NPV-1,1);
  iface  = 0;
  for ipar = 2:NPV;
    vorotmp = zeros(NVMX,2);
    itmp = 0;
    for ivo = 1:k;
      if(voroidx(ivo,ipar) == 1);
        itmp = itmp + 1;
        vorotmp(itmp,:) = [voro(ivo,1),voro(ivo,ipar)];
      end;
    end;
    niface(ipar-1) = itmp-1;
    if(niface(ipar-1) > 0);
      for ivo = 1:itmp-1;
        iface = iface + 1;
        ziface(iface) = vorotmp(ivo,1)+(vorotmp(ivo+1,1)-vorotmp(ivo,1))/2.;
        voroz(ivo,ipar-1,1) = ziface(iface);
        voroz(ivo,ipar-1,2) = vorotmp(ivo,2);
      end;
      voroz(itmp,ipar-1,2) = vorotmp(itmp,2);
    end;
    voroz(1,ipar-1,2) = vorotmp(1,2);
  end;

  %%
  %%  Find unique interfaces
  %%
  ntot = sum(niface);
  if(ntot > 0);
    ziface = sort(ziface);
    ziface = unique(ziface);
    nunique = length(ziface);
  else;
    nunique = 0;
  end;

if(nunique > 0);
  hiface = 0.;
  hiface(1) = ziface(1);
  for ivo = 2:nunique;
    hiface(ivo) = ziface(ivo)-ziface(ivo-1);
  end;

  partmp = 0.;
  partmp(1:nunique,1) = ziface(1:nunique);
  for ipar = 2:NPV;
    ivo = 1;
    if(niface(ipar-1) ~= 0);
      for ilay = 1:nunique;
        if(ziface(ilay) <= voroz(ivo,ipar-1,1));
          partmp(ilay,ipar) = voroz(ivo,ipar-1,2);
        else;
          ivo = ivo + 1;
          partmp(ilay,ipar) = voroz(ivo,ipar-1,2);
        end;
        if(ziface(ilay) >= voroz(niface(ipar-1),ipar-1,1));
          partmp(ilay+1:nunique+1,ipar) = voroz(ivo+1,ipar-1,2);
          break;
        end;
      end;
    else;
      partmp(1:nunique+1,ipar) = voroz(1,ipar-1,2);
    end;
  end;
  partmp(1:nunique,1) = ziface(1:nunique);
  par = 0.;
  for ilay = 1:nunique;
    par((ilay-1)*NPV+1:ilay*NPV) = partmp(ilay,:);
  end;
  par(nunique*NPV+1:((nunique+1)*NPV)-1) = partmp(nunique+1,2:end);
else;
  nunique = 0;
  par(1:NPV-1) = voro(1,2:end);
end;

return;
