#!/usr/bin/env python


import configparser
import imp
import json
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from pprint import pprint

import jsonschema
from jsonschema import Draft7Validator, validators


def import_config_default_smk():
    # inject VPIPE_BASEDIR into global namespace to make import below work:
    __builtins__.VPIPE_BASEDIR = "VPIPE_BASEDIR"

    # https://stackoverflow.com/questions/2601047/
    module_spec = spec_from_loader(
        "config_default", SourceFileLoader("config_default", "rules/config_default.smk")
    )
    config_default = module_from_spec(module_spec)
    module_spec.loader.exec_module(config_default)
    return config_default


def create_schema():

    config_default = import_config_default_smk()

    type_map = {
        str: "string",
        "str": "string",
        int: "integer",
        float: "number",
        bool: "boolean",
    }

    properties = {}

    try:
        config_default.VpipeConfig.__MEMBER_DEFAULT__
    except AttributeError:
        raise AttributeError(
            "looks like the config_default.smk was already rewritten"
        ) from None

    for section_name, section in config_default.VpipeConfig.__MEMBER_DEFAULT__.items():
        section_properties = {}
        for entry, record in section.items():
            section_properties[entry] = dict(
                type=type_map[record.type],
            )
            if record.value is not None:
                if section_name is "general" or entry != "threads":
                    section_properties[entry]["default"] = record.value
        properties[section_name] = dict(
            properties=section_properties,
            default={},
            type="object",
        )

    schema = {
        "type": "object",
        "$schema": "http://json-schema.org/draft-07/schema#",
        "properties": properties,
    }
    return schema


schema = create_schema()

with open("rules/config_schema.json", "w") as fh:
    json.dump(schema, fh)
