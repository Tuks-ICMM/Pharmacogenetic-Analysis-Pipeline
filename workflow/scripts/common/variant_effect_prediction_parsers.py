import logging
from typing import Union

from pandas import Series

def collectConsequenceTypes(apiResponse: dict) -> str:
    """
    Returns the most severe consequence for a variant.
    """
    if "most_severe_consequence" in apiResponse.keys():
        return apiResponse["most_severe_consequence"]
    else:
        return ""


def collectFeatureTypes(variant_effect_prediction: dict) -> str:
    """
    Returns the feature type recorded for a variant.
    """
    if "transcript_consequences" in variant_effect_prediction:
        return "Transcript"
    if "regulatory_feature_consequences" in variant_effect_prediction:
        return "Regulatory Feature"

def extractTranscriptConsequenceValue(vep_result: dict, field_name: str, log_instance_reference: str) -> Union[str,int]:
    """
    Returns the requested consequence type record for a variant.
    """
    _logger = logging.getLogger(log_instance_reference)
    # PROPERTY = list()
    _logger.info("Looking for %s field. Keys available to choose from: %s", field_name, ", ".join(vep_result.keys()))
    if "transcript_consequences" in vep_result:
        _logger.info("Transcript consequence identified. %s found.", len(vep_result["transcript_consequences"]))
        for consequence in vep_result["transcript_consequences"]:
            if field_name in consequence:
                _logger.info("%s field '%s' identified.", type(consequence), field_name)
                if isinstance(consequence[field_name], str):
                    _logger.info("Requested transcript consequence is a string.")
                    return consequence[field_name].replace("_", " ").capitalize()
                else: 
                    _logger.info("Requested transcript consequence is not a string.")
                    return consequence[field_name]
            else:
                return ""
    if "regulatory_feature_consequences" in vep_result:
        _logger.info("Regulatory feature consequence identified. %s found.", len(vep_result["regulatory_feature_consequences"]))
        for consequence in vep_result["regulatory_feature_consequences"]:
            if field_name in consequence:
                _logger.info("%s field '%s' identified.", type(consequence), field_name)
                if isinstance(consequence[field_name], str):
                    _logger.info("Requested regulatory feature consequence is a string.")
                    return consequence[field_name].replace("_", " ").capitalize()
                else: 
                    _logger.info("Requested regulatory feature consequence is not a string.")
                    return consequence[field_name]
            else:
                return ""
            
def extractMutFuncValues(vep_result: dict, field: str, calculation: str):
    if "transcript_consequences" in vep_result:
        for consequence in vep_result["transcript_consequences"]:
            if "mutfunc" in consequence and field in consequence["mutfunc"]:
                return vep_result["transcript_consequences"][0]["mutfunc"][field][calculation]
            else:
                return None

def extractUTRAnnotatorValues(vep_result: str, field: str):
    if "transcript_consequences" in vep_result:
        for consequence in vep_result["transcript_consequences"]:
            if field in consequence:
                return consequence[field]
            else:
                return None
    else:
        return None