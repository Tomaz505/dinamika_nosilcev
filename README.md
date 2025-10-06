# Dinamika nosilcev
Skripte za račun nelinearne dinamike nosilcev.

## Problemi
- Dinamika se računa za velike vrednosti mase
- Račun konvergira za velike časovne korake. Npr. $dt = 1.0 s$.
- Dinamika divergira za večje število KE.
- Premislim katere enote so najbolj primerje.
- Premisli številske vrednosti za $C$.
- Statika se izačuna ampak ob konstanti obtežbi alternira. (Tudi če je na začetku obremenjen element in nato ne več pomiki in hitrosti alternirajo predznake. V absolutnem so enaki.)
- Pri dinamiki robni pogoji sčasoma niso zagotovljeni.
- Deformacije imajo velike skokemed vrednostmi. $\kappa$ ima ustrezno obliko pod enakomernimi momenti če so končni elementi enake velikosti.

## Kontoliraj
- Kaj se zgodi z alterniranjem v statiki če uporabim za N nastavek, ki ne bi ohranjal energije? Poiskusi $R  =(R_1+R_2)/2$.
- Kako se račun obnaša z upoštevanjem $g$.

## Dopolnitve
- Dodaj točkovne obtežbe na robovih.
- Nekonzervativna obtežba (Potrebna linearizacija. Mislim, da samo rotacija v globalne komponente ki sledijo deformirani obliki elementa)
- Poišči bolši način računa uteži Lobatto integracije, če obstaja.
- Računski postopki za časovne končne elemente
- Sprostitve v krajiščih elementov.
- Pomična sprostitev v poljubni smeri.

## Spremembe
- Podatke pomikov in hitrosti spravi v $3d$ matriko. Mogoče celo samo hitrosti.
- Posplošitev določenih korakov v interacijskem računu, da bo omogočal različne metode.
- Premisli kako bi podajal obtežbo za večje število KE in, da za ničelno obtežbo ni potrebno spreminjati velikosti matrike.
