import os


class CsvGenerator:
    def __init__(self, lstring, sim_name: str, dirpath: str = 'outputs') -> None:
        self.lstring = lstring
        self.sim_name = sim_name
        self.dirpath = dirpath

    def metadata_to_csv(self, lsystem) -> None:
        out = open(os.path.join(self.dirpath, self.sim_name + '_metadata.csv'), 'w')
        out.write(f"GDD;{lsystem.derivationLength}\n")
        out.write(f"days;{lsystem.current_day}\n")
        out.write(f"plants;{lsystem.nb_plantes}\n")
        out.write(f"lines;{lsystem.NBlignes}\n")
        out.write(f"columns;{lsystem.NBcolonnes}\n")
        out.write(f"density;{lsystem.Densite}\n")
        out.write(f"spacing;{lsystem.Espacement}\n")
        out.write(f"start;{lsystem.sowing_date}\n")
        out.write(f"site;{lsystem.site}\n")
        out.write(f"tillers;{lsystem.option_tallage}\n")
        out.write(f"cutting_start;{lsystem.cutting_dates[0] if len(lsystem.cutting_dates) > 0 else None}\n")
        out.write(f"cutting_freq;{lsystem.cutting_dates[1] - lsystem.cutting_dates[0] if len(lsystem.cutting_dates) > 1 else None}\n")
        hcut = lsystem.ParamP[0]['HCOUP'] if len(set([i['HCOUP'] for i in lsystem.ParamP])) == 1 else 'variable'
        out.write(f"cutting_height;{hcut}\n")
        out.write(f"flowering;{lsystem.option_floraison}\n")
        out.write(f"senescence;{lsystem.option_senescence}\n")
        out.write(f"tiller_regression;{lsystem.option_tiller_regression}\n")
        out.write(f"carbon_regulation;{lsystem.option_morphogenetic_regulation_by_carbone}\n")
        out.close()

    def leaves_to_csv(self) -> None:
        out = open(os.path.join(self.dirpath, self.sim_name + '_feuilles.csv'), 'w')
        out.write("id_geno;id_plante;id_talle;id_rang;age;Agecroiss;Taillefeuille;Ymax;Taillelimbe;"
                  "Taillefinalelimbe;Taillegaine;Taillefinalegaine;Difftps;Phase;rapportK;coupe;Cutstatus;"
                  "angleinsert;angletal;surface_limbe;surface_gaine;biomass;Besoinencroiss;TailleEmergence;R\n")
        for mod in self.lstring:
            if mod.name in ('Feuille',):
                out.write(';'.join([str(mod[0].id_geno), str(mod[0].id_plante), str(mod[0].id_talle),
                                    str(mod[0].id_rang), str(mod[0].age), str(mod[0].Agecroiss),
                                    str(mod[0].Taillefeuille), str(mod[0].Ymax), str(mod[0].Taillelimbe),
                                    str(mod[0].Taillefinalelimbe), str(mod[0].Taillegaine),
                                    str(mod[0].Taillefinalegaine), str(mod[0].Difftps), str(mod[0].Phase),
                                    str(mod[0].rapportK), str(mod[0].coupe), str(mod[0].Cutstatus),
                                    str(mod[0].angleinsert), str(mod[0].angletal), str(mod[0].surface_limbe),
                                    str(mod[0].surface_gaine), str(mod[0].biomass), str(mod[0].Besoinencroiss),
                                    str(mod[0].TailleEmergence), str(mod[0].R)]) + '\n')
        out.close()

    def phytomers_to_csv(self) -> None:
        out = open(os.path.join(self.dirpath, self.sim_name + '_phytomeres.csv'), 'w')
        out.write("id_plante;id_talle;id_rang;age;angletal;organ;length\n")
        for mod in self.lstring:
            if mod.name in ('phytomere',):
                for rg, org in mod[0].organ_lengths.items():
                    for organ, length in org.items():
                        out.write(';'.join(
                            [str(mod[0].id_plante), str(mod[0].id_talle), str(rg), str(mod[0].age),
                             str(mod[0].angletal), str(organ), str(length)]) + '\n')
        out.close()

    def internodes_to_csv(self) -> None:
        out = open(os.path.join(self.dirpath, self.sim_name + '_entrenoeuds.csv'), 'w')
        out.write("id_plante;id_talle;id_rang;type;width;length;final_length\n")
        for mod in self.lstring:
            if mod.name in ('Entrenoeud',):
                out.write(';'.join([str(mod[0].id_plante), str(mod[0].id_talle), str(mod[0].id_rang),
                                    str(mod[0].internode_type), str(mod[0].width), str(mod[0].length),
                                    str(mod[0].final_length)]) + '\n')
        out.close()

    def apex_to_csv(self) -> None:
        out = open(os.path.join(self.dirpath, self.sim_name + '_apex.csv'), 'w')
        out.write("id_plante;id_talle;id_rang;id_geno;C;Premiecroiss;\n")
        for mod in self.lstring:
            if mod.name in ('apex',):
                out.write(';'.join([str(mod[0].id_plante), str(mod[0].id_talle), str(mod[0].id_rang),
                                    str(mod[0].id_geno), str(mod[0].C), str(mod[0].Premiecroiss)]) + '\n')
        out.close()
