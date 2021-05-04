import collections
import configparser
import os
import typing

import yaml  # TODO: record this dependency somewhere
import json
import jsonschema
from jsonschema import Draft7Validator, validators

__author__ = "Susana Posada-Cespedes"
__author__ = "David Seifert"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


# Class to parse config file
#
# the config object sets defaults for a number of parameters and
# validates parameters given by the user to be sensible


# https://python-jsonschema.readthedocs.io/en/latest/faq/#why-doesn-t-my-schema-s-default-property-set-the-default-on-my-instance

def extend_with_default(validator_class):
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(validator, properties, instance, schema):
        for property, subschema in properties.items():
            if "default" in subschema:
                default_value = subschema["default"]
                if isinstance(default_value, str):
                    default_value = default_value.format(**globals())
                instance.setdefault(property, default_value)

        not_allowed = set(instance) - set(properties)
        if not_allowed:
            message_ = "entry/ies {} not allowed.".format(
                ", ".join(sorted(not_allowed))
            )
            yield jsonschema.exceptions.ValidationError(message_)

        for error in validate_properties(
            validator,
            properties,
            instance,
            schema,
        ):
            yield error

    return validators.extend(
        validator_class,
        {"properties": set_defaults},
    )


DefaultValidatingDraft7Validator = extend_with_default(Draft7Validator)


VPIPE_CONFIG_OPTS = os.environ.get("VPIPE_CONFIG_OPTS") is not None


class _SectionWrapper:
    def __init__(self, config: "VpipeConfig", section):
        if not config.has_section(section):
            breakpoint()
            raise KeyError(
                "ERROR: Section '{}' is not a valid section!".format(section)
            )
        self._config = config
        self._section = section

    def __setitem__(self, option, value):
        self._config.set_option(self._section, option, value)

    def __getitem__(self, option):
        return self._config.get_option(self._section, option)


class VpipeConfig:
    "Class used to encapsulate the configuration properties used by V-pipe"

    def __init__(self, schema):
        # track all options (explicitly-set and defaults)
        self.__members = {}

        # track exclicitly-set options
        self._vpipe_configfile = configparser.ConfigParser()
        self._vpipe_configfile.read("vpipe.config")

        # TODO: rework whole config system (e.g. improve config merging)
        #with open(self.general["virus_config_file"]) as fd:
            #virus_config = yaml.safe_load(fd)

        #self.__members["virus_config"] = {}
        #for key, value in virus_config.items():
            #self.set_option("virus_config", key, value)
        config = {}
        for (name, section) in self._vpipe_configfile.items():
            config[name] = {}
            for (entry, value) in section.items():
                # fix different spellings for boolean values
                if value.lower() in ("true", "false"):
                    value = value.lower() == "true"
                config[name][entry] = value

        # DEFAULT entry is always created by configparser
        config.pop("DEFAULT", None)


        try:
            DefaultValidatingDraft7Validator(schema).validate(config)
        except jsonschema.exceptions.ValidationError as e:
            if len(e.path):
                section = e.path.pop()
                raise ValueError(f"error in section {section}: {e.message}")
            raise ValueError(f"error on top level: {e.message}")

        # inherit threads from general section if not specified by user:
        for (name, section) in config.items():
            if name == "general":
                continue
            if "threads" not in section:
                section["threads"] = config["general"]["threads"]

        self.__members = config

    def __getattr__(self, section_name) -> _SectionWrapper:
        return _SectionWrapper(self, section_name)

    def has_section(self, section):
        return section in self.__members

    def set_option(self, section, option, value):
        if section not in self.__members:
            raise KeyError(
                "ERROR: Section '{}' is not a valid section!".format(section)
            )
        self.__members[section][option] = value

        if section not in self._vpipe_configfile:
            self._vpipe_configfile[section] = {}
        # now add this explicitly set option
        self._vpipe_configfile[section][option] = value

    def get_option(self, section, option):
        if section not in self.__members:
            raise KeyError(
                "ERROR: Section '{}' is not a valid section!".format(section)
            )
        elif option not in self.__members[section]:
            raise ValueError(
                "ERROR: Section '{}' has no property '{}'!".format(section, option)
            )
        else:
            return self.__members[section][option]

    def write(self, outfile):
        self._vpipe_configfile.write(outfile)
