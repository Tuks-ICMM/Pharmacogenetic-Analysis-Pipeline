from json import dumps


class VariantFrequencies:
    """This object represents the database-derived variant frequencies for a given variant."""

    def __init__(
        self,
        af: float = 0.0,
        gnomade_nfe: float = 0.0,
        gnomadg_ami: float = 0.0,
        eas: float = 0.0,
        gnomadg_mid: float = 0.0,
        gnomade_eas: float = 0.0,
        gnomadg_oth: float = 0.0,
        gnomadg_eas: float = 0.0,
        gnomadg_amr: float = 0.0,
        gnomade_oth: float = 0.0,
        gnomade: float = 0.0,
        sas: float = 0.0,
        gnomadg: float = 0.0,
        gnomade_sas: float = 0.0,
        gnomadg_afr: float = 0.0,
        gnomadg_asj: float = 0.0,
        afr: float = 0.0,
        gnomade_afr: float = 0.0,
        eur: float = 0.0,
        gnomadg_nfe: float = 0.0,
        gnomadg_sas: float = 0.0,
        amr: float = 0.0,
        gnomade_fin: float = 0.0,
        gnomadg_fin: float = 0.0,
        gnomade_asj: float = 0.0,
        gnomade_amr: float = 0.0,
    ):
        self.af = af
        self.gnomade_nfe = gnomade_nfe
        self.gnomadg_ami = gnomadg_ami
        self.eas = eas
        self.gnomadg_mid = gnomadg_mid
        self.gnomade_eas = gnomade_eas
        self.gnomadg_oth = gnomadg_oth
        self.gnomadg_eas = gnomadg_eas
        self.gnomadg_amr = gnomadg_amr
        self.gnomade_oth = gnomade_oth
        self.gnomade = gnomade
        self.sas = sas
        self.gnomadg = gnomadg
        self.gnomade_sas = gnomade_sas
        self.gnomadg_afr = gnomadg_afr
        self.gnomadg_asj = gnomadg_asj
        self.afr = afr
        self.gnomade_afr = gnomade_afr
        self.eur = eur
        self.gnomadg_nfe = gnomadg_nfe
        self.gnomadg_sas = gnomadg_sas
        self.amr = amr
        self.gnomade_fin = gnomade_fin
        self.gnomadg_fin = gnomadg_fin
        self.gnomade_asj = gnomade_asj
        self.gnomade_amr = gnomade_amr

    def toDict(self):
        return self.__dict__
