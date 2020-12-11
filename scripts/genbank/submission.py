#!/usr/bin/env python
# mypy: ignore-errors
# -*- coding: utf-8 -*-

#
# Generated Thu Dec 10 13:22:26 2020 by generateDS.py version 2.37.8.
# Python 3.8.3 (default, Sep  7 2020, 14:20:01)  [GCC 7.5.0]
#
# Command line options:
#   ('-o', '/home/tautvydas/PycharmProjects/ps-bahrain-covid/scripts/genbank/genbank_submission.py')
#
# Command line arguments:
#   submission.xsd
#
# Command line:
#   /home/tautvydas/.pyenv/versions/3.8.3/envs/ps-bahrain-covid-38/bin/generateDS -o "/home/tautvydas/PycharmProjects/ps-bahrain-covid/scripts/genbank/genbank_submission.py" submission.xsd
#
# Current working directory (os.getcwd()):
#   ps-bahrain-covid
#

import sys

try:
    ModulenotfoundExp_ = ModuleNotFoundError
except NameError:
    ModulenotfoundExp_ = ImportError
from six.moves import zip_longest
import os
import re as re_
import base64
import datetime as datetime_
import decimal as decimal_

try:
    from lxml import etree as etree_
except ModulenotfoundExp_:
    from xml.etree import ElementTree as etree_


Validate_simpletypes_ = True
SaveElementTreeNode = True
if sys.version_info.major == 2:
    BaseStrType_ = basestring
else:
    BaseStrType_ = str


def parsexml_(infile, parser=None, **kwargs):
    if parser is None:
        # Use the lxml ElementTree compatible parser so that, e.g.,
        #   we ignore comments.
        try:
            parser = etree_.ETCompatXMLParser()
        except AttributeError:
            # fallback to xml.etree
            parser = etree_.XMLParser()
    try:
        if isinstance(infile, os.PathLike):
            infile = os.path.join(infile)
    except AttributeError:
        pass
    doc = etree_.parse(infile, parser=parser, **kwargs)
    return doc


def parsexmlstring_(instring, parser=None, **kwargs):
    if parser is None:
        # Use the lxml ElementTree compatible parser so that, e.g.,
        #   we ignore comments.
        try:
            parser = etree_.ETCompatXMLParser()
        except AttributeError:
            # fallback to xml.etree
            parser = etree_.XMLParser()
    element = etree_.fromstring(instring, parser=parser, **kwargs)
    return element


#
# Namespace prefix definition table (and other attributes, too)
#
# The module generatedsnamespaces, if it is importable, must contain
# a dictionary named GeneratedsNamespaceDefs.  This Python dictionary
# should map element type names (strings) to XML schema namespace prefix
# definitions.  The export method for any class for which there is
# a namespace prefix definition, will export that definition in the
# XML representation of that element.  See the export method of
# any generated element type class for an example of the use of this
# table.
# A sample table is:
#
#     # File: generatedsnamespaces.py
#
#     GenerateDSNamespaceDefs = {
#         "ElementtypeA": "http://www.xxx.com/namespaceA",
#         "ElementtypeB": "http://www.xxx.com/namespaceB",
#     }
#
# Additionally, the generatedsnamespaces module can contain a python
# dictionary named GenerateDSNamespaceTypePrefixes that associates element types
# with the namespace prefixes that are to be added to the
# "xsi:type" attribute value.  See the exportAttributes method of
# any generated element type and the generation of "xsi:type" for an
# example of the use of this table.
# An example table:
#
#     # File: generatedsnamespaces.py
#
#     GenerateDSNamespaceTypePrefixes = {
#         "ElementtypeC": "aaa:",
#         "ElementtypeD": "bbb:",
#     }
#

try:
    from generatedsnamespaces import GenerateDSNamespaceDefs as GenerateDSNamespaceDefs_
except ModulenotfoundExp_:
    GenerateDSNamespaceDefs_ = {}
try:
    from generatedsnamespaces import GenerateDSNamespaceTypePrefixes as GenerateDSNamespaceTypePrefixes_
except ModulenotfoundExp_:
    GenerateDSNamespaceTypePrefixes_ = {}

#
# You can replace the following class definition by defining an
# importable module named "generatedscollector" containing a class
# named "GdsCollector".  See the default class definition below for
# clues about the possible content of that class.
#
try:
    from generatedscollector import GdsCollector as GdsCollector_
except ModulenotfoundExp_:

    class GdsCollector_(object):
        def __init__(self, messages=None):
            if messages is None:
                self.messages = []
            else:
                self.messages = messages

        def add_message(self, msg):
            self.messages.append(msg)

        def get_messages(self):
            return self.messages

        def clear_messages(self):
            self.messages = []

        def print_messages(self):
            for msg in self.messages:
                print("Warning: {}".format(msg))

        def write_messages(self, outstream):
            for msg in self.messages:
                outstream.write("Warning: {}\n".format(msg))


#
# The super-class for enum types
#

try:
    from enum import Enum
except ModulenotfoundExp_:
    Enum = object

#
# The root super-class for element type classes
#
# Calls to the methods in these classes are generated by generateDS.py.
# You can replace these methods by re-implementing the following class
#   in a module named generatedssuper.py.

try:
    from generatedssuper import GeneratedsSuper
except ModulenotfoundExp_ as exp:

    class GeneratedsSuper(object):
        __hash__ = object.__hash__
        tzoff_pattern = re_.compile(r"(\+|-)((0\d|1[0-3]):[0-5]\d|14:00)$")

        class _FixedOffsetTZ(datetime_.tzinfo):
            def __init__(self, offset, name):
                self.__offset = datetime_.timedelta(minutes=offset)
                self.__name = name

            def utcoffset(self, dt):
                return self.__offset

            def tzname(self, dt):
                return self.__name

            def dst(self, dt):
                return None

        def gds_format_string(self, input_data, input_name=""):
            return input_data

        def gds_parse_string(self, input_data, node=None, input_name=""):
            return input_data

        def gds_validate_string(self, input_data, node=None, input_name=""):
            if not input_data:
                return ""
            else:
                return input_data

        def gds_format_base64(self, input_data, input_name=""):
            return base64.b64encode(input_data)

        def gds_validate_base64(self, input_data, node=None, input_name=""):
            return input_data

        def gds_format_integer(self, input_data, input_name=""):
            return "%d" % input_data

        def gds_parse_integer(self, input_data, node=None, input_name=""):
            try:
                ival = int(input_data)
            except (TypeError, ValueError) as exp:
                raise_parse_error(node, "Requires integer value: %s" % exp)
            return ival

        def gds_validate_integer(self, input_data, node=None, input_name=""):
            try:
                value = int(input_data)
            except (TypeError, ValueError):
                raise_parse_error(node, "Requires integer value")
            return value

        def gds_format_integer_list(self, input_data, input_name=""):
            if len(input_data) > 0 and not isinstance(input_data[0], BaseStrType_):
                input_data = [str(s) for s in input_data]
            return "%s" % " ".join(input_data)

        def gds_validate_integer_list(self, input_data, node=None, input_name=""):
            values = input_data.split()
            for value in values:
                try:
                    int(value)
                except (TypeError, ValueError):
                    raise_parse_error(node, "Requires sequence of integer valuess")
            return values

        def gds_format_float(self, input_data, input_name=""):
            return ("%.15f" % input_data).rstrip("0")

        def gds_parse_float(self, input_data, node=None, input_name=""):
            try:
                fval_ = float(input_data)
            except (TypeError, ValueError) as exp:
                raise_parse_error(node, "Requires float or double value: %s" % exp)
            return fval_

        def gds_validate_float(self, input_data, node=None, input_name=""):
            try:
                value = float(input_data)
            except (TypeError, ValueError):
                raise_parse_error(node, "Requires float value")
            return value

        def gds_format_float_list(self, input_data, input_name=""):
            if len(input_data) > 0 and not isinstance(input_data[0], BaseStrType_):
                input_data = [str(s) for s in input_data]
            return "%s" % " ".join(input_data)

        def gds_validate_float_list(self, input_data, node=None, input_name=""):
            values = input_data.split()
            for value in values:
                try:
                    float(value)
                except (TypeError, ValueError):
                    raise_parse_error(node, "Requires sequence of float values")
            return values

        def gds_format_decimal(self, input_data, input_name=""):
            return_value = "%s" % input_data
            if "." in return_value:
                return_value = return_value.rstrip("0")
                if return_value.endswith("."):
                    return_value = return_value.rstrip(".")
            return return_value

        def gds_parse_decimal(self, input_data, node=None, input_name=""):
            try:
                decimal_value = decimal_.Decimal(input_data)
            except (TypeError, ValueError):
                raise_parse_error(node, "Requires decimal value")
            return decimal_value

        def gds_validate_decimal(self, input_data, node=None, input_name=""):
            try:
                value = decimal_.Decimal(input_data)
            except (TypeError, ValueError):
                raise_parse_error(node, "Requires decimal value")
            return value

        def gds_format_decimal_list(self, input_data, input_name=""):
            if len(input_data) > 0 and not isinstance(input_data[0], BaseStrType_):
                input_data = [str(s) for s in input_data]
            return " ".join([self.gds_format_decimal(item) for item in input_data])

        def gds_validate_decimal_list(self, input_data, node=None, input_name=""):
            values = input_data.split()
            for value in values:
                try:
                    decimal_.Decimal(value)
                except (TypeError, ValueError):
                    raise_parse_error(node, "Requires sequence of decimal values")
            return values

        def gds_format_double(self, input_data, input_name=""):
            return "%s" % input_data

        def gds_parse_double(self, input_data, node=None, input_name=""):
            try:
                fval_ = float(input_data)
            except (TypeError, ValueError) as exp:
                raise_parse_error(node, "Requires double or float value: %s" % exp)
            return fval_

        def gds_validate_double(self, input_data, node=None, input_name=""):
            try:
                value = float(input_data)
            except (TypeError, ValueError):
                raise_parse_error(node, "Requires double or float value")
            return value

        def gds_format_double_list(self, input_data, input_name=""):
            if len(input_data) > 0 and not isinstance(input_data[0], BaseStrType_):
                input_data = [str(s) for s in input_data]
            return "%s" % " ".join(input_data)

        def gds_validate_double_list(self, input_data, node=None, input_name=""):
            values = input_data.split()
            for value in values:
                try:
                    float(value)
                except (TypeError, ValueError):
                    raise_parse_error(node, "Requires sequence of double or float values")
            return values

        def gds_format_boolean(self, input_data, input_name=""):
            return ("%s" % input_data).lower()

        def gds_parse_boolean(self, input_data, node=None, input_name=""):
            if input_data in ("true", "1"):
                bval = True
            elif input_data in ("false", "0"):
                bval = False
            else:
                raise_parse_error(node, "Requires boolean value")
            return bval

        def gds_validate_boolean(self, input_data, node=None, input_name=""):
            if input_data not in (
                True,
                1,
                False,
                0,
            ):
                raise_parse_error(node, "Requires boolean value " "(one of True, 1, False, 0)")
            return input_data

        def gds_format_boolean_list(self, input_data, input_name=""):
            if len(input_data) > 0 and not isinstance(input_data[0], BaseStrType_):
                input_data = [str(s) for s in input_data]
            return "%s" % " ".join(input_data)

        def gds_validate_boolean_list(self, input_data, node=None, input_name=""):
            values = input_data.split()
            for value in values:
                if value not in (
                    True,
                    1,
                    False,
                    0,
                ):
                    raise_parse_error(node, "Requires sequence of boolean values " "(one of True, 1, False, 0)")
            return values

        def gds_validate_datetime(self, input_data, node=None, input_name=""):
            return input_data

        def gds_format_datetime(self, input_data, input_name=""):
            if input_data.microsecond == 0:
                _svalue = "%04d-%02d-%02dT%02d:%02d:%02d" % (
                    input_data.year,
                    input_data.month,
                    input_data.day,
                    input_data.hour,
                    input_data.minute,
                    input_data.second,
                )
            else:
                _svalue = "%04d-%02d-%02dT%02d:%02d:%02d.%s" % (
                    input_data.year,
                    input_data.month,
                    input_data.day,
                    input_data.hour,
                    input_data.minute,
                    input_data.second,
                    ("%f" % (float(input_data.microsecond) / 1000000))[2:],
                )
            if input_data.tzinfo is not None:
                tzoff = input_data.tzinfo.utcoffset(input_data)
                if tzoff is not None:
                    total_seconds = tzoff.seconds + (86400 * tzoff.days)
                    if total_seconds == 0:
                        _svalue += "Z"
                    else:
                        if total_seconds < 0:
                            _svalue += "-"
                            total_seconds *= -1
                        else:
                            _svalue += "+"
                        hours = total_seconds // 3600
                        minutes = (total_seconds - (hours * 3600)) // 60
                        _svalue += "{0:02d}:{1:02d}".format(hours, minutes)
            return _svalue

        @classmethod
        def gds_parse_datetime(cls, input_data):
            tz = None
            if input_data[-1] == "Z":
                tz = GeneratedsSuper._FixedOffsetTZ(0, "UTC")
                input_data = input_data[:-1]
            else:
                results = GeneratedsSuper.tzoff_pattern.search(input_data)
                if results is not None:
                    tzoff_parts = results.group(2).split(":")
                    tzoff = int(tzoff_parts[0]) * 60 + int(tzoff_parts[1])
                    if results.group(1) == "-":
                        tzoff *= -1
                    tz = GeneratedsSuper._FixedOffsetTZ(tzoff, results.group(0))
                    input_data = input_data[:-6]
            time_parts = input_data.split(".")
            if len(time_parts) > 1:
                micro_seconds = int(float("0." + time_parts[1]) * 1000000)
                input_data = "%s.%s" % (
                    time_parts[0],
                    "{}".format(micro_seconds).rjust(6, "0"),
                )
                dt = datetime_.datetime.strptime(input_data, "%Y-%m-%dT%H:%M:%S.%f")
            else:
                dt = datetime_.datetime.strptime(input_data, "%Y-%m-%dT%H:%M:%S")
            dt = dt.replace(tzinfo=tz)
            return dt

        def gds_validate_date(self, input_data, node=None, input_name=""):
            return input_data

        def gds_format_date(self, input_data, input_name=""):
            _svalue = "%04d-%02d-%02d" % (
                input_data.year,
                input_data.month,
                input_data.day,
            )
            try:
                if input_data.tzinfo is not None:
                    tzoff = input_data.tzinfo.utcoffset(input_data)
                    if tzoff is not None:
                        total_seconds = tzoff.seconds + (86400 * tzoff.days)
                        if total_seconds == 0:
                            _svalue += "Z"
                        else:
                            if total_seconds < 0:
                                _svalue += "-"
                                total_seconds *= -1
                            else:
                                _svalue += "+"
                            hours = total_seconds // 3600
                            minutes = (total_seconds - (hours * 3600)) // 60
                            _svalue += "{0:02d}:{1:02d}".format(hours, minutes)
            except AttributeError:
                pass
            return _svalue

        @classmethod
        def gds_parse_date(cls, input_data):
            tz = None
            if input_data[-1] == "Z":
                tz = GeneratedsSuper._FixedOffsetTZ(0, "UTC")
                input_data = input_data[:-1]
            else:
                results = GeneratedsSuper.tzoff_pattern.search(input_data)
                if results is not None:
                    tzoff_parts = results.group(2).split(":")
                    tzoff = int(tzoff_parts[0]) * 60 + int(tzoff_parts[1])
                    if results.group(1) == "-":
                        tzoff *= -1
                    tz = GeneratedsSuper._FixedOffsetTZ(tzoff, results.group(0))
                    input_data = input_data[:-6]
            dt = datetime_.datetime.strptime(input_data, "%Y-%m-%d")
            dt = dt.replace(tzinfo=tz)
            return dt.date()

        def gds_validate_time(self, input_data, node=None, input_name=""):
            return input_data

        def gds_format_time(self, input_data, input_name=""):
            if input_data.microsecond == 0:
                _svalue = "%02d:%02d:%02d" % (
                    input_data.hour,
                    input_data.minute,
                    input_data.second,
                )
            else:
                _svalue = "%02d:%02d:%02d.%s" % (
                    input_data.hour,
                    input_data.minute,
                    input_data.second,
                    ("%f" % (float(input_data.microsecond) / 1000000))[2:],
                )
            if input_data.tzinfo is not None:
                tzoff = input_data.tzinfo.utcoffset(input_data)
                if tzoff is not None:
                    total_seconds = tzoff.seconds + (86400 * tzoff.days)
                    if total_seconds == 0:
                        _svalue += "Z"
                    else:
                        if total_seconds < 0:
                            _svalue += "-"
                            total_seconds *= -1
                        else:
                            _svalue += "+"
                        hours = total_seconds // 3600
                        minutes = (total_seconds - (hours * 3600)) // 60
                        _svalue += "{0:02d}:{1:02d}".format(hours, minutes)
            return _svalue

        def gds_validate_simple_patterns(self, patterns, target):
            # pat is a list of lists of strings/patterns.
            # The target value must match at least one of the patterns
            # in order for the test to succeed.
            found1 = True
            for patterns1 in patterns:
                found2 = False
                for patterns2 in patterns1:
                    mo = re_.search(patterns2, target)
                    if mo is not None and len(mo.group(0)) == len(target):
                        found2 = True
                        break
                if not found2:
                    found1 = False
                    break
            return found1

        @classmethod
        def gds_parse_time(cls, input_data):
            tz = None
            if input_data[-1] == "Z":
                tz = GeneratedsSuper._FixedOffsetTZ(0, "UTC")
                input_data = input_data[:-1]
            else:
                results = GeneratedsSuper.tzoff_pattern.search(input_data)
                if results is not None:
                    tzoff_parts = results.group(2).split(":")
                    tzoff = int(tzoff_parts[0]) * 60 + int(tzoff_parts[1])
                    if results.group(1) == "-":
                        tzoff *= -1
                    tz = GeneratedsSuper._FixedOffsetTZ(tzoff, results.group(0))
                    input_data = input_data[:-6]
            if len(input_data.split(".")) > 1:
                dt = datetime_.datetime.strptime(input_data, "%H:%M:%S.%f")
            else:
                dt = datetime_.datetime.strptime(input_data, "%H:%M:%S")
            dt = dt.replace(tzinfo=tz)
            return dt.time()

        def gds_check_cardinality_(self, value, input_name, min_occurs=0, max_occurs=1, required=None):
            if value is None:
                length = 0
            elif isinstance(value, list):
                length = len(value)
            else:
                length = 1
            if required is not None:
                if required and length < 1:
                    self.gds_collector_.add_message(
                        "Required value {}{} is missing".format(input_name, self.gds_get_node_lineno_())
                    )
            if length < min_occurs:
                self.gds_collector_.add_message(
                    "Number of values for {}{} is below "
                    "the minimum allowed, "
                    "expected at least {}, found {}".format(input_name, self.gds_get_node_lineno_(), min_occurs, length)
                )
            elif length > max_occurs:
                self.gds_collector_.add_message(
                    "Number of values for {}{} is above "
                    "the maximum allowed, "
                    "expected at most {}, found {}".format(input_name, self.gds_get_node_lineno_(), max_occurs, length)
                )

        def gds_validate_builtin_ST_(
            self, validator, value, input_name, min_occurs=None, max_occurs=None, required=None
        ):
            if value is not None:
                try:
                    validator(value, input_name=input_name)
                except GDSParseError as parse_error:
                    self.gds_collector_.add_message(str(parse_error))

        def gds_validate_defined_ST_(
            self, validator, value, input_name, min_occurs=None, max_occurs=None, required=None
        ):
            if value is not None:
                try:
                    validator(value)
                except GDSParseError as parse_error:
                    self.gds_collector_.add_message(str(parse_error))

        def gds_str_lower(self, instring):
            return instring.lower()

        def get_path_(self, node):
            path_list = []
            self.get_path_list_(node, path_list)
            path_list.reverse()
            path = "/".join(path_list)
            return path

        Tag_strip_pattern_ = re_.compile(r"\{.*\}")

        def get_path_list_(self, node, path_list):
            if node is None:
                return
            tag = GeneratedsSuper.Tag_strip_pattern_.sub("", node.tag)
            if tag:
                path_list.append(tag)
            self.get_path_list_(node.getparent(), path_list)

        def get_class_obj_(self, node, default_class=None):
            class_obj1 = default_class
            if "xsi" in node.nsmap:
                classname = node.get("{%s}type" % node.nsmap["xsi"])
                if classname is not None:
                    names = classname.split(":")
                    if len(names) == 2:
                        classname = names[1]
                    class_obj2 = globals().get(classname)
                    if class_obj2 is not None:
                        class_obj1 = class_obj2
            return class_obj1

        def gds_build_any(self, node, type_name=None):
            # provide default value in case option --disable-xml is used.
            content = ""
            content = etree_.tostring(node, encoding="unicode")
            return content

        @classmethod
        def gds_reverse_node_mapping(cls, mapping):
            return dict(((v, k) for k, v in mapping.items()))

        @staticmethod
        def gds_encode(instring):
            if sys.version_info.major == 2:
                if ExternalEncoding:
                    encoding = ExternalEncoding
                else:
                    encoding = "utf-8"
                return instring.encode(encoding)
            else:
                return instring

        @staticmethod
        def convert_unicode(instring):
            if isinstance(instring, str):
                result = quote_xml(instring)
            elif sys.version_info.major == 2 and isinstance(instring, unicode):
                result = quote_xml(instring).encode("utf8")
            else:
                result = GeneratedsSuper.gds_encode(str(instring))
            return result

        def __eq__(self, other):
            def excl_select_objs_(obj):
                return obj[0] != "parent_object_" and obj[0] != "gds_collector_"

            if type(self) != type(other):
                return False
            return all(
                x == y
                for x, y in zip_longest(
                    filter(excl_select_objs_, self.__dict__.items()), filter(excl_select_objs_, other.__dict__.items())
                )
            )

        def __ne__(self, other):
            return not self.__eq__(other)

        # Django ETL transform hooks.
        def gds_djo_etl_transform(self):
            pass

        def gds_djo_etl_transform_db_obj(self, dbobj):
            pass

        # SQLAlchemy ETL transform hooks.
        def gds_sqa_etl_transform(self):
            return 0, None

        def gds_sqa_etl_transform_db_obj(self, dbobj):
            pass

        def gds_get_node_lineno_(self):
            if hasattr(self, "gds_elementtree_node_") and self.gds_elementtree_node_ is not None:
                return " near line {}".format(self.gds_elementtree_node_.sourceline)
            else:
                return ""

    def getSubclassFromModule_(module, class_):
        """Get the subclass of a class from a specific module."""
        name = class_.__name__ + "Sub"
        if hasattr(module, name):
            return getattr(module, name)
        else:
            return None


#
# If you have installed IPython you can uncomment and use the following.
# IPython is available from http://ipython.scipy.org/.
#

## from IPython.Shell import IPShellEmbed
## args = ''
## ipshell = IPShellEmbed(args,
##     banner = 'Dropping into IPython',
##     exit_msg = 'Leaving Interpreter, back to program.')

# Then use the following line where and when you want to drop into the
# IPython shell:
#    ipshell('<some message> -- Entering ipshell.\nHit Ctrl-D to exit')

#
# Globals
#

ExternalEncoding = ""
# Set this to false in order to deactivate during export, the use of
# name space prefixes captured from the input document.
UseCapturedNS_ = True
CapturedNsmap_ = {}
Tag_pattern_ = re_.compile(r"({.*})?(.*)")
String_cleanup_pat_ = re_.compile(r"[\n\r\s]+")
Namespace_extract_pat_ = re_.compile(r"{(.*)}(.*)")
CDATA_pattern_ = re_.compile(r"<!\[CDATA\[.*?\]\]>", re_.DOTALL)

# Change this to redirect the generated superclass module to use a
# specific subclass module.
CurrentSubclassModule_ = None

#
# Support/utility functions.
#


def showIndent(outfile, level, pretty_print=True):
    if pretty_print:
        for idx in range(level):
            outfile.write("    ")


def quote_xml(inStr):
    "Escape markup chars, but do not modify CDATA sections."
    if not inStr:
        return ""
    s1 = isinstance(inStr, BaseStrType_) and inStr or "%s" % inStr
    s2 = ""
    pos = 0
    matchobjects = CDATA_pattern_.finditer(s1)
    for mo in matchobjects:
        s3 = s1[pos : mo.start()]
        s2 += quote_xml_aux(s3)
        s2 += s1[mo.start() : mo.end()]
        pos = mo.end()
    s3 = s1[pos:]
    s2 += quote_xml_aux(s3)
    return s2


def quote_xml_aux(inStr):
    s1 = inStr.replace("&", "&amp;")
    s1 = s1.replace("<", "&lt;")
    s1 = s1.replace(">", "&gt;")
    return s1


def quote_attrib(inStr):
    s1 = isinstance(inStr, BaseStrType_) and inStr or "%s" % inStr
    s1 = s1.replace("&", "&amp;")
    s1 = s1.replace("<", "&lt;")
    s1 = s1.replace(">", "&gt;")
    if '"' in s1:
        if "'" in s1:
            s1 = '"%s"' % s1.replace('"', "&quot;")
        else:
            s1 = "'%s'" % s1
    else:
        s1 = '"%s"' % s1
    return s1


def quote_python(inStr):
    s1 = inStr
    if s1.find("'") == -1:
        if s1.find("\n") == -1:
            return "'%s'" % s1
        else:
            return "'''%s'''" % s1
    else:
        if s1.find('"') != -1:
            s1 = s1.replace('"', '\\"')
        if s1.find("\n") == -1:
            return '"%s"' % s1
        else:
            return '"""%s"""' % s1


def get_all_text_(node):
    if node.text is not None:
        text = node.text
    else:
        text = ""
    for child in node:
        if child.tail is not None:
            text += child.tail
    return text


def find_attr_value_(attr_name, node):
    attrs = node.attrib
    attr_parts = attr_name.split(":")
    value = None
    if len(attr_parts) == 1:
        value = attrs.get(attr_name)
    elif len(attr_parts) == 2:
        prefix, name = attr_parts
        namespace = node.nsmap.get(prefix)
        if namespace is not None:
            value = attrs.get(
                "{%s}%s"
                % (
                    namespace,
                    name,
                )
            )
    return value


def encode_str_2_3(instr):
    return instr


class GDSParseError(Exception):
    pass


def raise_parse_error(node, msg):
    if node is not None:
        msg = "%s (element %s/line %d)" % (
            msg,
            node.tag,
            node.sourceline,
        )
    raise GDSParseError(msg)


class MixedContainer:
    # Constants for category:
    CategoryNone = 0
    CategoryText = 1
    CategorySimple = 2
    CategoryComplex = 3
    # Constants for content_type:
    TypeNone = 0
    TypeText = 1
    TypeString = 2
    TypeInteger = 3
    TypeFloat = 4
    TypeDecimal = 5
    TypeDouble = 6
    TypeBoolean = 7
    TypeBase64 = 8

    def __init__(self, category, content_type, name, value):
        self.category = category
        self.content_type = content_type
        self.name = name
        self.value = value

    def getCategory(self):
        return self.category

    def getContenttype(self, content_type):
        return self.content_type

    def getValue(self):
        return self.value

    def getName(self):
        return self.name

    def export(self, outfile, level, name, namespace, pretty_print=True):
        if self.category == MixedContainer.CategoryText:
            # Prevent exporting empty content as empty lines.
            if self.value.strip():
                outfile.write(self.value)
        elif self.category == MixedContainer.CategorySimple:
            self.exportSimple(outfile, level, name)
        else:  # category == MixedContainer.CategoryComplex
            self.value.export(outfile, level, namespace, name_=name, pretty_print=pretty_print)

    def exportSimple(self, outfile, level, name):
        if self.content_type == MixedContainer.TypeString:
            outfile.write("<%s>%s</%s>" % (self.name, self.value, self.name))
        elif self.content_type == MixedContainer.TypeInteger or self.content_type == MixedContainer.TypeBoolean:
            outfile.write("<%s>%d</%s>" % (self.name, self.value, self.name))
        elif self.content_type == MixedContainer.TypeFloat or self.content_type == MixedContainer.TypeDecimal:
            outfile.write("<%s>%f</%s>" % (self.name, self.value, self.name))
        elif self.content_type == MixedContainer.TypeDouble:
            outfile.write("<%s>%g</%s>" % (self.name, self.value, self.name))
        elif self.content_type == MixedContainer.TypeBase64:
            outfile.write("<%s>%s</%s>" % (self.name, base64.b64encode(self.value), self.name))

    def to_etree(self, element, mapping_=None, nsmap_=None):
        if self.category == MixedContainer.CategoryText:
            # Prevent exporting empty content as empty lines.
            if self.value.strip():
                if len(element) > 0:
                    if element[-1].tail is None:
                        element[-1].tail = self.value
                    else:
                        element[-1].tail += self.value
                else:
                    if element.text is None:
                        element.text = self.value
                    else:
                        element.text += self.value
        elif self.category == MixedContainer.CategorySimple:
            subelement = etree_.SubElement(element, "%s" % self.name)
            subelement.text = self.to_etree_simple()
        else:  # category == MixedContainer.CategoryComplex
            self.value.to_etree(element)

    def to_etree_simple(self, mapping_=None, nsmap_=None):
        if self.content_type == MixedContainer.TypeString:
            text = self.value
        elif self.content_type == MixedContainer.TypeInteger or self.content_type == MixedContainer.TypeBoolean:
            text = "%d" % self.value
        elif self.content_type == MixedContainer.TypeFloat or self.content_type == MixedContainer.TypeDecimal:
            text = "%f" % self.value
        elif self.content_type == MixedContainer.TypeDouble:
            text = "%g" % self.value
        elif self.content_type == MixedContainer.TypeBase64:
            text = "%s" % base64.b64encode(self.value)
        return text

    def exportLiteral(self, outfile, level, name):
        if self.category == MixedContainer.CategoryText:
            showIndent(outfile, level)
            outfile.write(
                'model_.MixedContainer(%d, %d, "%s", "%s"),\n'
                % (self.category, self.content_type, self.name, self.value)
            )
        elif self.category == MixedContainer.CategorySimple:
            showIndent(outfile, level)
            outfile.write(
                'model_.MixedContainer(%d, %d, "%s", "%s"),\n'
                % (self.category, self.content_type, self.name, self.value)
            )
        else:  # category == MixedContainer.CategoryComplex
            showIndent(outfile, level)
            outfile.write(
                'model_.MixedContainer(%d, %d, "%s",\n'
                % (
                    self.category,
                    self.content_type,
                    self.name,
                )
            )
            self.value.exportLiteral(outfile, level + 1)
            showIndent(outfile, level)
            outfile.write(")\n")


class MemberSpec_(object):
    def __init__(self, name="", data_type="", container=0, optional=0, child_attrs=None, choice=None):
        self.name = name
        self.data_type = data_type
        self.container = container
        self.child_attrs = child_attrs
        self.choice = choice
        self.optional = optional

    def set_name(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def set_data_type(self, data_type):
        self.data_type = data_type

    def get_data_type_chain(self):
        return self.data_type

    def get_data_type(self):
        if isinstance(self.data_type, list):
            if len(self.data_type) > 0:
                return self.data_type[-1]
            else:
                return "xs:string"
        else:
            return self.data_type

    def set_container(self, container):
        self.container = container

    def get_container(self):
        return self.container

    def set_child_attrs(self, child_attrs):
        self.child_attrs = child_attrs

    def get_child_attrs(self):
        return self.child_attrs

    def set_choice(self, choice):
        self.choice = choice

    def get_choice(self):
        return self.choice

    def set_optional(self, optional):
        self.optional = optional

    def get_optional(self):
        return self.optional


def _cast(typ, value):
    if typ is None or value is None:
        return value
    return typ(value)


#
# Data representation classes.
#


class DataTypeType(str, Enum):
    """DataType is what we can process:
    submitter-xml,project-core-xml,biosample-xml,genbank-FF,genbank-seqsubmit,
    sra-experiment.xml,sra-run.xml, etc
    full list must be specified in submission database. Branded with the
    version."""

    AUTODETECTXML = "autodetect-xml"  # XML is autodected by top-level element and XML schema reference
    GENERICDATA = "generic-data"  # This is a data file, the type of which will be determined by destination archive
    PHENOTYPETABLE = "phenotype-table"  # Phenotype table
    SRASTUDYXMLV_1 = "sra-study-xml-v1"
    SRAEXPERIMENTXMLV_1 = "sra-experiment-xml-v1"
    SRASAMPLEXMLV_1 = "sra-sample-xml-v1"
    SRARUNXMLV_1 = "sra-run-xml-v1"
    SRAANALYSISXMLV_1 = "sra-analysis-xml-v1"
    SRASTUDYXMLV_2 = "sra-study-xml-v2"
    SRAEXPERIMENTXMLV_2 = "sra-experiment-xml-v2"
    SRASAMPLEXMLV_2 = "sra-sample-xml-v2"
    SRARUNXMLV_2 = "sra-run-xml-v2"
    SRAANALYSISXMLV_2 = "sra-analysis-xml-v2"
    SRARUN_454_NATIVE = "sra-run-454_native"
    SRARUNBAM = "sra-run-bam"
    SRARUN_COMPLETE_GENOMICS_NATIVE = "sra-run-CompleteGenomics_native"
    SRARUNFASTQ = "sra-run-fastq"
    SRARUN_HELICOS_NATIVE = "sra-run-Helicos_native"
    SRARUN_PAC_BIO_HDF_5 = "sra-run-PacBio_HDF5"
    SRARUNSFF = "sra-run-sff"
    SRARUNSO_LI_D_NATIVE = "sra-run-SOLiD_native"
    SRARUNSRF = "sra-run-srf"
    PROJECTCOREXMLV_1 = "project-core-xml-v1"
    WGSCONTIGSSQN = "wgs-contigs-sqn"
    WGSUNPLACEDSCAFFOLDSAGP = "wgs-unplaced-scaffolds-agp"
    WGSCONTIGREPLICONDESCR = "wgs-contig-replicon-descr"
    WGSAGPREPLICONDESCR = "wgs-agp-replicon-descr"
    WGSLOCCHRTOREPLICON = "wgs-loc-chr-to-replicon"
    WGSREPLICONFROMCONTIGSAGP = "wgs-replicon-from-contigs-agp"
    WGSSCAFFOLDFROMCONTIGSAGP = "wgs-scaffold-from-contigs-agp"
    WGSREPLICONFROMSCAFFOLDSAGP = "wgs-replicon-from-scaffolds-agp"
    WGSUNLOCALIZEDSCAFFOLDSAGP = "wgs-unlocalized-scaffolds-agp"
    WGSUNLOCSCAFFOLDTOREPLICON = "wgs-unloc-scaffold-to-replicon"
    WGSASSEMBLYSQN = "wgs-assembly-sqn"
    WGSASSEMBLYFASTA = "wgs-assembly-fasta"
    WGSCONTIGSFASTA = "wgs-contigs-fasta"
    WGSAGP = "wgs-agp"
    WGSPLACEMENT = "wgs-placement"
    ENAWGSFLATFILE = "ena-wgs-flatfile"
    DDBJWGSFLATFILE = "ddbj-wgs-flatfile"
    WGSFLATFILEPREPROCESSREPORT = "wgs-flatfile-preprocess-report"
    TSASEQSUBMITSQN = "tsa-seqsubmit-sqn"
    COMPLETEGENOMESANNOTATEDSQN = "complete-genomes-annotated-sqn"
    COMPLETEGENOMESANNOTATESQN = "complete-genomes-annotate-sqn"
    COMPLETEGENOMESANNOTATEFASTA = "complete-genomes-annotate-fasta"
    COMPLETEGENOMESANNOTATETEMPLATE = "complete-genomes-annotate-template"
    COMPLETEGENOMESREPLICON = "complete-genomes-replicon"
    GENBANKSQN = "genbank-sqn"
    GENBANKSUBMISSIONPACKAGE = "genbank-submission-package"
    GENBANKBARCODETAR = "genbank-barcode-tar"
    GENBANKSEQUENCESFASTA = "genbank-sequences-fasta"
    GENBANKSRCMODSTBL = "genbank-srcmods-tbl"
    GENBANKNCBILINKTBL = "genbank-ncbi-link-tbl"
    GENBANKSEQUENCESFILTEREDFASTA = "genbank-sequences-filtered-fasta"
    GENBANKSEQUENCESFASTAVALXML = "genbank-sequences-fastaval-xml"
    GENBANKSRCMODSFILTEREDTBL = "genbank-srcmods-filtered-tbl"
    GENBANKSEQUENCESREPORTTXT = "genbank-sequences-report-txt"
    GENBANKSEQUENCESREPORTTBL = "genbank-sequences-report-tbl"
    GENBANKTOOLSVERSIONSXML = "genbank-tools-versions-xml"
    GENBANKFEATURESTABLE = "genbank-features-table"
    GENBANKFEATURESFILTEREDTABLE = "genbank-features-filtered-table"
    METHYLATIONDATA = "methylation-data"
    SEQUENCESFASTA = "sequences-fasta"
    BIONANOCMAP = "bionano-cmap"
    BIONANOCOORD = "bionano-coord"
    BIONANOXMAP = "bionano-xmap"
    BIONANOSMAP = "bionano-smap"
    BIONANOBNX = "bionano-bnx"
    SEQUIN = "sequin"
    ANTIBIOGRAM = "antibiogram"
    MODIFICATIONSCSV = "modifications.csv"
    MODIFICATIONSGFF = "modifications.gff"
    MOTIFSGFF = "motifs.gff"
    MOTIF_SUMMARYCSV = "motif_summary.csv"
    BIOSAMPLETBLV_2_0 = "biosample-tbl-v2.0"
    ANTIBIOGRAMTBLV_1_0 = "antibiogram-tbl-v1.0"


class DbTypeType(str, Enum):
    E_PMC = "ePMC"
    E_PUBMED = "ePubmed"
    E_DOI = "eDOI"
    E_NOT_AVAILABLE = "eNotAvailable"  # Use for Unpublished, in press or other databases.


class content_encodingType(str, Enum):
    """How data is encoded (or how to decode it) E.g. - plain or base64"""

    PLAIN = "plain"  # Plain text
    BASE_64 = "base64"  # Base64-encoded binary


class roleType(str, Enum):
    """Role of the ogranization in submission - owner of the data or just a
    participant. It is expected that there is one owner of the submission
    data."""

    OWNER = "owner"
    PARTICIPANT = "participant"


class statusType(str, Enum):
    E_PUBLISHED = "ePublished"
    E_IN_PRESS = "eInPress"
    E_UNPUBLISHED = "eUnpublished"


class typeArchive(str, Enum):
    """One of INSDC archives."""

    NCBI = "NCBI"
    EBI = "EBI"
    DDBJ = "DDBJ"


class typeTargetDb(str, Enum):
    """Supported target databases"""

    BIO_PROJECT = "BioProject"
    BIO_SAMPLE = "BioSample"
    CLINVAR = "clinvar"  # ClinVar Submission
    DB_GA_P = "dbGaP"  # dbGaP phenotypic data
    WGS = "WGS"  # WGS GenBank submission
    VARIATION = "variation"  # Variation Organization submission
    VARIATION_SUBMISSION = "variation_submission"  # Variation File submission
    GTR = "GTR"  # GTR Lab submission
    TSA = "TSA"  # TSA GenBank submission
    COMPLETE_GENOMES = "CompleteGenomes"  # CompleteGenomes GenBank submission
    DB_VAR = "dbVar"
    SRA = "SRA"
    SRAEXPERIMENT = "SRA.experiment"
    SRARUN = "SRA.run"
    SP = "SP"
    PGAP = "PGAP"  # PGAP annotations in Submission Portal
    GEN_BANK = "GenBank"  # GenBank Web submissions processor
    SUP_FILES = "SupFiles"
    EST = "EST"  # EST GenBank submission
    GSS = "GSS"  # GSS GenBank submission
    TPA = "TPA"  # TPA GenBank submission


class typeType(str, Enum):
    """Organization type : center, institute, consortium or medical lab"""

    CONSORTIUM = "consortium"
    CENTER = "center"
    INSTITUTE = "institute"
    LAB = "lab"


class Submission(GeneratedsSuper):
    """Reference to the existing submission in the case of re-submit."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        schema_version=None,
        resubmit_of=None,
        submitted=None,
        last_update=None,
        status=None,
        submission_id=None,
        Description=None,
        Action=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.schema_version = _cast(None, schema_version)
        self.schema_version_nsprefix_ = None
        self.resubmit_of = _cast(None, resubmit_of)
        self.resubmit_of_nsprefix_ = None
        if isinstance(submitted, BaseStrType_):
            initvalue_ = datetime_.datetime.strptime(submitted, "%Y-%m-%d").date()
        else:
            initvalue_ = submitted
        self.submitted = initvalue_
        if isinstance(last_update, BaseStrType_):
            initvalue_ = datetime_.datetime.strptime(last_update, "%Y-%m-%d").date()
        else:
            initvalue_ = last_update
        self.last_update = initvalue_
        self.status = _cast(None, status)
        self.status_nsprefix_ = None
        self.submission_id = _cast(None, submission_id)
        self.submission_id_nsprefix_ = None
        self.Description = Description
        self.Description_nsprefix_ = None
        if Action is None:
            self.Action = []
        else:
            self.Action = Action
        self.Action_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, Submission)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if Submission.subclass:
            return Submission.subclass(*args_, **kwargs_)
        else:
            return Submission(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Description(self):
        return self.Description

    def set_Description(self, Description):
        self.Description = Description

    def get_Action(self):
        return self.Action

    def set_Action(self, Action):
        self.Action = Action

    def add_Action(self, value):
        self.Action.append(value)

    def insert_Action_at(self, index, value):
        self.Action.insert(index, value)

    def replace_Action_at(self, index, value):
        self.Action[index] = value

    def get_schema_version(self):
        return self.schema_version

    def set_schema_version(self, schema_version):
        self.schema_version = schema_version

    def get_resubmit_of(self):
        return self.resubmit_of

    def set_resubmit_of(self, resubmit_of):
        self.resubmit_of = resubmit_of

    def get_submitted(self):
        return self.submitted

    def set_submitted(self, submitted):
        self.submitted = submitted

    def get_last_update(self):
        return self.last_update

    def set_last_update(self, last_update):
        self.last_update = last_update

    def get_status(self):
        return self.status

    def set_status(self, status):
        self.status = status

    def get_submission_id(self):
        return self.submission_id

    def set_submission_id(self, submission_id):
        self.submission_id = submission_id

    def hasContent_(self):
        if self.Description is not None or self.Action:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="Submission", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("Submission")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "Submission":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="Submission")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="Submission", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="Submission"):
        if self.schema_version is not None and "schema_version" not in already_processed:
            already_processed.add("schema_version")
            outfile.write(
                " schema_version=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(quote_attrib(self.schema_version), input_name="schema_version")
                    ),
                )
            )
        if self.resubmit_of is not None and "resubmit_of" not in already_processed:
            already_processed.add("resubmit_of")
            outfile.write(
                " resubmit_of=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.resubmit_of), input_name="resubmit_of")),)
            )
        if self.submitted is not None and "submitted" not in already_processed:
            already_processed.add("submitted")
            outfile.write(' submitted="%s"' % self.gds_format_date(self.submitted, input_name="submitted"))
        if self.last_update is not None and "last_update" not in already_processed:
            already_processed.add("last_update")
            outfile.write(' last_update="%s"' % self.gds_format_date(self.last_update, input_name="last_update"))
        if self.status is not None and "status" not in already_processed:
            already_processed.add("status")
            outfile.write(
                " status=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.status), input_name="status")),)
            )
        if self.submission_id is not None and "submission_id" not in already_processed:
            already_processed.add("submission_id")
            outfile.write(
                " submission_id=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(quote_attrib(self.submission_id), input_name="submission_id")
                    ),
                )
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="Submission",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Description is not None:
            namespaceprefix_ = (
                self.Description_nsprefix_ + ":" if (UseCapturedNS_ and self.Description_nsprefix_) else ""
            )
            self.Description.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Description", pretty_print=pretty_print
            )
        for Action_ in self.Action:
            namespaceprefix_ = self.Action_nsprefix_ + ":" if (UseCapturedNS_ and self.Action_nsprefix_) else ""
            Action_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Action", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("schema_version", node)
        if value is not None and "schema_version" not in already_processed:
            already_processed.add("schema_version")
            self.schema_version = value
        value = find_attr_value_("resubmit_of", node)
        if value is not None and "resubmit_of" not in already_processed:
            already_processed.add("resubmit_of")
            self.resubmit_of = value
        value = find_attr_value_("submitted", node)
        if value is not None and "submitted" not in already_processed:
            already_processed.add("submitted")
            try:
                self.submitted = self.gds_parse_date(value)
            except ValueError as exp:
                raise ValueError("Bad date attribute (submitted): %s" % exp)
        value = find_attr_value_("last_update", node)
        if value is not None and "last_update" not in already_processed:
            already_processed.add("last_update")
            try:
                self.last_update = self.gds_parse_date(value)
            except ValueError as exp:
                raise ValueError("Bad date attribute (last_update): %s" % exp)
        value = find_attr_value_("status", node)
        if value is not None and "status" not in already_processed:
            already_processed.add("status")
            self.status = value
        value = find_attr_value_("submission_id", node)
        if value is not None and "submission_id" not in already_processed:
            already_processed.add("submission_id")
            self.submission_id = value
            self.submission_id = " ".join(self.submission_id.split())

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Description":
            obj_ = DescriptionType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Description = obj_
            obj_.original_tagname_ = "Description"
        elif nodeName_ == "Action":
            obj_ = ActionType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Action.append(obj_)
            obj_.original_tagname_ = "Action"


# end class Submission


class typeFile(GeneratedsSuper):
    """Path to the file relative from the location of submission XML. Note - at
    least one of (i) file_path, (ii) file_id or (iii) md5 needs to be
    present to correctly address the file.
    FileTrack file id - unique and more preferred way to address the file.
    Works for files already in FileTrack.
    FileTrack cloud file id - unique id for file in the cloud, which metadata
    is stored in FileTrack.
    Purpose of md5 is two-fold: to verify file content or to address the file.
    Using of crc32 is infer compare to md5 and expected to be used in internal
    processing only.
    Standart content type - e.g. text/xml, etc."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        file_path=None,
        file_id=None,
        cloud_id=None,
        md5=None,
        crc32=None,
        content_type=None,
        DataType=None,
        extensiontype_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.file_path = _cast(None, file_path)
        self.file_path_nsprefix_ = None
        self.file_id = _cast(None, file_id)
        self.file_id_nsprefix_ = None
        self.cloud_id = _cast(None, cloud_id)
        self.cloud_id_nsprefix_ = None
        self.md5 = _cast(None, md5)
        self.md5_nsprefix_ = None
        self.crc32 = _cast(None, crc32)
        self.crc32_nsprefix_ = None
        self.content_type = _cast(None, content_type)
        self.content_type_nsprefix_ = None
        self.DataType = DataType
        if self.DataType is not None:
            self.validate_DataTypeType(self.DataType)
        self.DataType_nsprefix_ = None
        self.extensiontype_ = extensiontype_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeFile)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeFile.subclass:
            return typeFile.subclass(*args_, **kwargs_)
        else:
            return typeFile(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_DataType(self):
        return self.DataType

    def set_DataType(self, DataType):
        self.DataType = DataType

    def get_file_path(self):
        return self.file_path

    def set_file_path(self, file_path):
        self.file_path = file_path

    def get_file_id(self):
        return self.file_id

    def set_file_id(self, file_id):
        self.file_id = file_id

    def get_cloud_id(self):
        return self.cloud_id

    def set_cloud_id(self, cloud_id):
        self.cloud_id = cloud_id

    def get_md5(self):
        return self.md5

    def set_md5(self, md5):
        self.md5 = md5

    def get_crc32(self):
        return self.crc32

    def set_crc32(self, crc32):
        self.crc32 = crc32

    def get_content_type(self):
        return self.content_type

    def set_content_type(self, content_type):
        self.content_type = content_type

    def get_extensiontype_(self):
        return self.extensiontype_

    def set_extensiontype_(self, extensiontype_):
        self.extensiontype_ = extensiontype_

    def validate_DataTypeType(self, value):
        result = True
        # Validate type DataTypeType, a restriction on xs:string.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, str):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (str)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            value = value
            enumerations = [
                "autodetect-xml",
                "generic-data",
                "phenotype-table",
                "sra-study-xml-v1",
                "sra-experiment-xml-v1",
                "sra-sample-xml-v1",
                "sra-run-xml-v1",
                "sra-analysis-xml-v1",
                "sra-study-xml-v2",
                "sra-experiment-xml-v2",
                "sra-sample-xml-v2",
                "sra-run-xml-v2",
                "sra-analysis-xml-v2",
                "sra-run-454_native",
                "sra-run-bam",
                "sra-run-CompleteGenomics_native",
                "sra-run-fastq",
                "sra-run-Helicos_native",
                "sra-run-PacBio_HDF5",
                "sra-run-sff",
                "sra-run-SOLiD_native",
                "sra-run-srf",
                "project-core-xml-v1",
                "wgs-contigs-sqn",
                "wgs-unplaced-scaffolds-agp",
                "wgs-contig-replicon-descr",
                "wgs-agp-replicon-descr",
                "wgs-loc-chr-to-replicon",
                "wgs-replicon-from-contigs-agp",
                "wgs-scaffold-from-contigs-agp",
                "wgs-replicon-from-scaffolds-agp",
                "wgs-unlocalized-scaffolds-agp",
                "wgs-unloc-scaffold-to-replicon",
                "wgs-assembly-sqn",
                "wgs-assembly-fasta",
                "wgs-contigs-fasta",
                "wgs-agp",
                "wgs-placement",
                "ena-wgs-flatfile",
                "ddbj-wgs-flatfile",
                "wgs-flatfile-preprocess-report",
                "tsa-seqsubmit-sqn",
                "complete-genomes-annotated-sqn",
                "complete-genomes-annotate-sqn",
                "complete-genomes-annotate-fasta",
                "complete-genomes-annotate-template",
                "complete-genomes-replicon",
                "genbank-sqn",
                "genbank-submission-package",
                "genbank-barcode-tar",
                "genbank-sequences-fasta",
                "genbank-srcmods-tbl",
                "genbank-ncbi-link-tbl",
                "genbank-sequences-filtered-fasta",
                "genbank-sequences-fastaval-xml",
                "genbank-srcmods-filtered-tbl",
                "genbank-sequences-report-txt",
                "genbank-sequences-report-tbl",
                "genbank-tools-versions-xml",
                "genbank-features-table",
                "genbank-features-filtered-table",
                "methylation-data",
                "sequences-fasta",
                "bionano-cmap",
                "bionano-coord",
                "bionano-xmap",
                "bionano-smap",
                "bionano-bnx",
                "sequin",
                "antibiogram",
                "modifications.csv",
                "modifications.gff",
                "motifs.gff",
                "motif_summary.csv",
                "biosample-tbl-v2.0",
                "antibiogram-tbl-v1.0",
            ]
            if value not in enumerations:
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s does not match xsd enumeration restriction on DataTypeType'
                    % {"value": encode_str_2_3(value), "lineno": lineno}
                )
                result = False
        return result

    def hasContent_(self):
        if self.DataType is not None:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="typeFile", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeFile")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeFile":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeFile")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeFile", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="typeFile"):
        if self.file_path is not None and "file_path" not in already_processed:
            already_processed.add("file_path")
            outfile.write(
                " file_path=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.file_path), input_name="file_path")),)
            )
        if self.file_id is not None and "file_id" not in already_processed:
            already_processed.add("file_id")
            outfile.write(
                " file_id=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.file_id), input_name="file_id")),)
            )
        if self.cloud_id is not None and "cloud_id" not in already_processed:
            already_processed.add("cloud_id")
            outfile.write(
                " cloud_id=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.cloud_id), input_name="cloud_id")),)
            )
        if self.md5 is not None and "md5" not in already_processed:
            already_processed.add("md5")
            outfile.write(
                " md5=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.md5), input_name="md5")),)
            )
        if self.crc32 is not None and "crc32" not in already_processed:
            already_processed.add("crc32")
            outfile.write(
                " crc32=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.crc32), input_name="crc32")),)
            )
        if self.content_type is not None and "content_type" not in already_processed:
            already_processed.add("content_type")
            outfile.write(
                " content_type=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.content_type), input_name="content_type")),)
            )
        if self.extensiontype_ is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            outfile.write(' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"')
            if ":" not in self.extensiontype_:
                imported_ns_type_prefix_ = GenerateDSNamespaceTypePrefixes_.get(self.extensiontype_, "")
                outfile.write(' xsi:type="%s%s"' % (imported_ns_type_prefix_, self.extensiontype_))
            else:
                outfile.write(' xsi:type="%s"' % self.extensiontype_)

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="typeFile",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.DataType is not None:
            namespaceprefix_ = self.DataType_nsprefix_ + ":" if (UseCapturedNS_ and self.DataType_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sDataType>%s</%sDataType>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.DataType), input_name="DataType")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("file_path", node)
        if value is not None and "file_path" not in already_processed:
            already_processed.add("file_path")
            self.file_path = value
        value = find_attr_value_("file_id", node)
        if value is not None and "file_id" not in already_processed:
            already_processed.add("file_id")
            self.file_id = value
        value = find_attr_value_("cloud_id", node)
        if value is not None and "cloud_id" not in already_processed:
            already_processed.add("cloud_id")
            self.cloud_id = value
        value = find_attr_value_("md5", node)
        if value is not None and "md5" not in already_processed:
            already_processed.add("md5")
            self.md5 = value
        value = find_attr_value_("crc32", node)
        if value is not None and "crc32" not in already_processed:
            already_processed.add("crc32")
            self.crc32 = value
        value = find_attr_value_("content_type", node)
        if value is not None and "content_type" not in already_processed:
            already_processed.add("content_type")
            self.content_type = value
        value = find_attr_value_("xsi:type", node)
        if value is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            self.extensiontype_ = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "DataType":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "DataType")
            value_ = self.gds_validate_string(value_, node, "DataType")
            self.DataType = value_
            self.DataType_nsprefix_ = child_.prefix
            # validate type DataTypeType
            self.validate_DataTypeType(self.DataType)


# end class typeFile


class typeInlineData(GeneratedsSuper):
    """This is inline data to be embedded into the submissionEither XML or
    base64/plain text data)Optional name of the data objectData model of
    the data objectContent type - what is it - XML, text, binary, etcHow
    data is encoded (or how to decode it) E.g. - plain or base64"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        name=None,
        data_model=None,
        content_type=None,
        content_encoding=None,
        XmlContent=None,
        DataContent=None,
        extensiontype_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.name = _cast(None, name)
        self.name_nsprefix_ = None
        self.data_model = _cast(None, data_model)
        self.data_model_nsprefix_ = None
        self.content_type = _cast(None, content_type)
        self.content_type_nsprefix_ = None
        self.content_encoding = _cast(None, content_encoding)
        self.content_encoding_nsprefix_ = None
        self.XmlContent = XmlContent
        self.XmlContent_nsprefix_ = None
        self.DataContent = DataContent
        self.DataContent_nsprefix_ = None
        self.extensiontype_ = extensiontype_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeInlineData)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeInlineData.subclass:
            return typeInlineData.subclass(*args_, **kwargs_)
        else:
            return typeInlineData(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_XmlContent(self):
        return self.XmlContent

    def set_XmlContent(self, XmlContent):
        self.XmlContent = XmlContent

    def get_DataContent(self):
        return self.DataContent

    def set_DataContent(self, DataContent):
        self.DataContent = DataContent

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

    def get_data_model(self):
        return self.data_model

    def set_data_model(self, data_model):
        self.data_model = data_model

    def get_content_type(self):
        return self.content_type

    def set_content_type(self, content_type):
        self.content_type = content_type

    def get_content_encoding(self):
        return self.content_encoding

    def set_content_encoding(self, content_encoding):
        self.content_encoding = content_encoding

    def get_extensiontype_(self):
        return self.extensiontype_

    def set_extensiontype_(self, extensiontype_):
        self.extensiontype_ = extensiontype_

    def validate_content_encodingType(self, value):
        # Validate type content_encodingType, a restriction on xs:string.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, str):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (str)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            value = value
            enumerations = ["plain", "base64"]
            if value not in enumerations:
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s does not match xsd enumeration restriction on content_encodingType'
                    % {"value": encode_str_2_3(value), "lineno": lineno}
                )
                result = False

    def hasContent_(self):
        if self.XmlContent is not None or self.DataContent is not None:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="typeInlineData", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeInlineData")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeInlineData":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeInlineData")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeInlineData", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="typeInlineData"):
        if self.name is not None and "name" not in already_processed:
            already_processed.add("name")
            outfile.write(
                " name=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.name), input_name="name")),)
            )
        if self.data_model is not None and "data_model" not in already_processed:
            already_processed.add("data_model")
            outfile.write(
                " data_model=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.data_model), input_name="data_model")),)
            )
        if self.content_type is not None and "content_type" not in already_processed:
            already_processed.add("content_type")
            outfile.write(
                " content_type=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.content_type), input_name="content_type")),)
            )
        if self.content_encoding is not None and "content_encoding" not in already_processed:
            already_processed.add("content_encoding")
            outfile.write(
                " content_encoding=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(quote_attrib(self.content_encoding), input_name="content_encoding")
                    ),
                )
            )
        if self.extensiontype_ is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            outfile.write(' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"')
            if ":" not in self.extensiontype_:
                imported_ns_type_prefix_ = GenerateDSNamespaceTypePrefixes_.get(self.extensiontype_, "")
                outfile.write(' xsi:type="%s%s"' % (imported_ns_type_prefix_, self.extensiontype_))
            else:
                outfile.write(' xsi:type="%s"' % self.extensiontype_)

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="typeInlineData",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.XmlContent is not None:
            namespaceprefix_ = self.XmlContent_nsprefix_ + ":" if (UseCapturedNS_ and self.XmlContent_nsprefix_) else ""
            self.XmlContent.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="XmlContent", pretty_print=pretty_print
            )
        if self.DataContent is not None:
            namespaceprefix_ = (
                self.DataContent_nsprefix_ + ":" if (UseCapturedNS_ and self.DataContent_nsprefix_) else ""
            )
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sDataContent>%s</%sDataContent>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.DataContent), input_name="DataContent")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("name", node)
        if value is not None and "name" not in already_processed:
            already_processed.add("name")
            self.name = value
        value = find_attr_value_("data_model", node)
        if value is not None and "data_model" not in already_processed:
            already_processed.add("data_model")
            self.data_model = value
        value = find_attr_value_("content_type", node)
        if value is not None and "content_type" not in already_processed:
            already_processed.add("content_type")
            self.content_type = value
        value = find_attr_value_("content_encoding", node)
        if value is not None and "content_encoding" not in already_processed:
            already_processed.add("content_encoding")
            self.content_encoding = value
            self.validate_content_encodingType(self.content_encoding)  # validate type content_encodingType
        value = find_attr_value_("xsi:type", node)
        if value is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            self.extensiontype_ = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "XmlContent":
            obj_ = XmlContentType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.XmlContent = obj_
            obj_.original_tagname_ = "XmlContent"
        elif nodeName_ == "DataContent":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "DataContent")
            value_ = self.gds_validate_string(value_, node, "DataContent")
            self.DataContent = value_
            self.DataContent_nsprefix_ = child_.prefix


# end class typeInlineData


class typeAccount(GeneratedsSuper):
    """Submission Portal account iddeprecateddeprecated"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, account_id=None, user_name=None, authority=None, Contact=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.account_id = _cast(None, account_id)
        self.account_id_nsprefix_ = None
        self.user_name = _cast(None, user_name)
        self.user_name_nsprefix_ = None
        self.authority = _cast(None, authority)
        self.authority_nsprefix_ = None
        self.Contact = Contact
        self.Contact_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeAccount)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeAccount.subclass:
            return typeAccount.subclass(*args_, **kwargs_)
        else:
            return typeAccount(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Contact(self):
        return self.Contact

    def set_Contact(self, Contact):
        self.Contact = Contact

    def get_account_id(self):
        return self.account_id

    def set_account_id(self, account_id):
        self.account_id = account_id

    def get_user_name(self):
        return self.user_name

    def set_user_name(self, user_name):
        self.user_name = user_name

    def get_authority(self):
        return self.authority

    def set_authority(self, authority):
        self.authority = authority

    def hasContent_(self):
        if self.Contact is not None:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeAccount",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeAccount")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeAccount":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeAccount")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeAccount", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="typeAccount"):
        if self.account_id is not None and "account_id" not in already_processed:
            already_processed.add("account_id")
            outfile.write(
                " account_id=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.account_id), input_name="account_id")),)
            )
        if self.user_name is not None and "user_name" not in already_processed:
            already_processed.add("user_name")
            outfile.write(
                " user_name=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.user_name), input_name="user_name")),)
            )
        if self.authority is not None and "authority" not in already_processed:
            already_processed.add("authority")
            outfile.write(
                " authority=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.authority), input_name="authority")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeAccount",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Contact is not None:
            namespaceprefix_ = self.Contact_nsprefix_ + ":" if (UseCapturedNS_ and self.Contact_nsprefix_) else ""
            self.Contact.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Contact", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("account_id", node)
        if value is not None and "account_id" not in already_processed:
            already_processed.add("account_id")
            self.account_id = value
        value = find_attr_value_("user_name", node)
        if value is not None and "user_name" not in already_processed:
            already_processed.add("user_name")
            self.user_name = value
        value = find_attr_value_("authority", node)
        if value is not None and "authority" not in already_processed:
            already_processed.add("authority")
            self.authority = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Contact":
            obj_ = typeContactInfo.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Contact = obj_
            obj_.original_tagname_ = "Contact"


# end class typeAccount


class typeOrganization(GeneratedsSuper):
    """Organization for the submissionOrganization type : center, institute,
    consortium or medical lab
    Role of the ogranization in submission - owner of the data or just a
    participant. It is expected that there is one owner of the submission
    data.
    In case we want to exchange organization listURL of the organization
    website."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        type_=None,
        role=None,
        org_id=None,
        url=None,
        group_id=None,
        Name=None,
        Address=None,
        Contact=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.type_ = _cast(None, type_)
        self.type__nsprefix_ = None
        self.role = _cast(None, role)
        self.role_nsprefix_ = None
        self.org_id = _cast(int, org_id)
        self.org_id_nsprefix_ = None
        self.url = _cast(None, url)
        self.url_nsprefix_ = None
        self.group_id = _cast(None, group_id)
        self.group_id_nsprefix_ = None
        self.Name = Name
        self.Name_nsprefix_ = None
        self.Address = Address
        self.Address_nsprefix_ = None
        if Contact is None:
            self.Contact = []
        else:
            self.Contact = Contact
        self.Contact_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeOrganization)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeOrganization.subclass:
            return typeOrganization.subclass(*args_, **kwargs_)
        else:
            return typeOrganization(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Name(self):
        return self.Name

    def set_Name(self, Name):
        self.Name = Name

    def get_Address(self):
        return self.Address

    def set_Address(self, Address):
        self.Address = Address

    def get_Contact(self):
        return self.Contact

    def set_Contact(self, Contact):
        self.Contact = Contact

    def add_Contact(self, value):
        self.Contact.append(value)

    def insert_Contact_at(self, index, value):
        self.Contact.insert(index, value)

    def replace_Contact_at(self, index, value):
        self.Contact[index] = value

    def get_type(self):
        return self.type_

    def set_type(self, type_):
        self.type_ = type_

    def get_role(self):
        return self.role

    def set_role(self, role):
        self.role = role

    def get_org_id(self):
        return self.org_id

    def set_org_id(self, org_id):
        self.org_id = org_id

    def get_url(self):
        return self.url

    def set_url(self, url):
        self.url = url

    def get_group_id(self):
        return self.group_id

    def set_group_id(self, group_id):
        self.group_id = group_id

    def validate_typeType(self, value):
        # Validate type typeType, a restriction on xs:string.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, str):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (str)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            value = value
            enumerations = ["consortium", "center", "institute", "lab"]
            if value not in enumerations:
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s does not match xsd enumeration restriction on typeType'
                    % {"value": encode_str_2_3(value), "lineno": lineno}
                )
                result = False

    def validate_roleType(self, value):
        # Validate type roleType, a restriction on xs:string.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, str):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (str)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            value = value
            enumerations = ["owner", "participant"]
            if value not in enumerations:
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s does not match xsd enumeration restriction on roleType'
                    % {"value": encode_str_2_3(value), "lineno": lineno}
                )
                result = False

    def hasContent_(self):
        if self.Name is not None or self.Address is not None or self.Contact:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeOrganization",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeOrganization")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeOrganization":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeOrganization")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeOrganization", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="typeOrganization"):
        if self.type_ is not None and "type_" not in already_processed:
            already_processed.add("type_")
            outfile.write(
                " type=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.type_), input_name="type")),)
            )
        if self.role is not None and "role" not in already_processed:
            already_processed.add("role")
            outfile.write(
                " role=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.role), input_name="role")),)
            )
        if self.org_id is not None and "org_id" not in already_processed:
            already_processed.add("org_id")
            outfile.write(' org_id="%s"' % self.gds_format_integer(self.org_id, input_name="org_id"))
        if self.url is not None and "url" not in already_processed:
            already_processed.add("url")
            outfile.write(
                " url=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.url), input_name="url")),)
            )
        if self.group_id is not None and "group_id" not in already_processed:
            already_processed.add("group_id")
            outfile.write(
                " group_id=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.group_id), input_name="group_id")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeOrganization",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Name is not None:
            namespaceprefix_ = self.Name_nsprefix_ + ":" if (UseCapturedNS_ and self.Name_nsprefix_) else ""
            self.Name.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Name", pretty_print=pretty_print
            )
        if self.Address is not None:
            namespaceprefix_ = self.Address_nsprefix_ + ":" if (UseCapturedNS_ and self.Address_nsprefix_) else ""
            self.Address.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Address", pretty_print=pretty_print
            )
        for Contact_ in self.Contact:
            namespaceprefix_ = self.Contact_nsprefix_ + ":" if (UseCapturedNS_ and self.Contact_nsprefix_) else ""
            Contact_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Contact", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("type", node)
        if value is not None and "type" not in already_processed:
            already_processed.add("type")
            self.type_ = value
            self.validate_typeType(self.type_)  # validate type typeType
        value = find_attr_value_("role", node)
        if value is not None and "role" not in already_processed:
            already_processed.add("role")
            self.role = value
            self.validate_roleType(self.role)  # validate type roleType
        value = find_attr_value_("org_id", node)
        if value is not None and "org_id" not in already_processed:
            already_processed.add("org_id")
            self.org_id = self.gds_parse_integer(value, node, "org_id")
            if self.org_id <= 0:
                raise_parse_error(node, "Invalid PositiveInteger")
        value = find_attr_value_("url", node)
        if value is not None and "url" not in already_processed:
            already_processed.add("url")
            self.url = value
        value = find_attr_value_("group_id", node)
        if value is not None and "group_id" not in already_processed:
            already_processed.add("group_id")
            self.group_id = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Name":
            obj_ = NameType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Name = obj_
            obj_.original_tagname_ = "Name"
        elif nodeName_ == "Address":
            obj_ = typeAddress.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Address = obj_
            obj_.original_tagname_ = "Address"
        elif nodeName_ == "Contact":
            obj_ = typeContactInfo.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Contact.append(obj_)
            obj_.original_tagname_ = "Contact"


# end class typeOrganization


class typeFileAttribute(GeneratedsSuper):
    """Named attributes, attached to the file. This way submitter can attribute
    file to BioProject, BioSample, etc."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, name=None, valueOf_=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.name = _cast(None, name)
        self.name_nsprefix_ = None
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeFileAttribute)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeFileAttribute.subclass:
            return typeFileAttribute.subclass(*args_, **kwargs_)
        else:
            return typeFileAttribute(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def hasContent_(self):
        if 1 if type(self.valueOf_) in [int, float] else self.valueOf_:
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="", namespacedef_="", name_="typeFileAttribute", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeFileAttribute")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeFileAttribute":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeFileAttribute")
        if self.hasContent_():
            outfile.write(">")
            outfile.write(self.convert_unicode(self.valueOf_))
            self.exportChildren(
                outfile,
                level + 1,
                namespaceprefix_,
                namespacedef_,
                name_="typeFileAttribute",
                pretty_print=pretty_print,
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="typeFileAttribute"):
        if self.name is not None and "name" not in already_processed:
            already_processed.add("name")
            outfile.write(
                " name=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.name), input_name="name")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="typeFileAttribute",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("name", node)
        if value is not None and "name" not in already_processed:
            already_processed.add("name")
            self.name = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class typeFileAttribute


class typeFileAttributeRefId(GeneratedsSuper):
    """Named attributes, attached to the file. This way submitter can attribute
    file to BioProject, BioSample, etc."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, name=None, RefId=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.name = _cast(None, name)
        self.name_nsprefix_ = None
        self.RefId = RefId
        self.RefId_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeFileAttributeRefId)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeFileAttributeRefId.subclass:
            return typeFileAttributeRefId.subclass(*args_, **kwargs_)
        else:
            return typeFileAttributeRefId(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_RefId(self):
        return self.RefId

    def set_RefId(self, RefId):
        self.RefId = RefId

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

    def hasContent_(self):
        if self.RefId is not None:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeFileAttributeRefId",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeFileAttributeRefId")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeFileAttributeRefId":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeFileAttributeRefId")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile,
                level + 1,
                namespaceprefix_,
                namespacedef_,
                name_="typeFileAttributeRefId",
                pretty_print=pretty_print,
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="typeFileAttributeRefId"):
        if self.name is not None and "name" not in already_processed:
            already_processed.add("name")
            outfile.write(
                " name=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.name), input_name="name")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeFileAttributeRefId",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.RefId is not None:
            namespaceprefix_ = self.RefId_nsprefix_ + ":" if (UseCapturedNS_ and self.RefId_nsprefix_) else ""
            self.RefId.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="RefId", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("name", node)
        if value is not None and "name" not in already_processed:
            already_processed.add("name")
            self.name = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "RefId":
            obj_ = typeRefId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.RefId = obj_
            obj_.original_tagname_ = "RefId"


# end class typeFileAttributeRefId


class typeSequenceData(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, Sequence=None, AuthorSet=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if Sequence is None:
            self.Sequence = []
        else:
            self.Sequence = Sequence
        self.Sequence_nsprefix_ = None
        self.AuthorSet = AuthorSet
        self.AuthorSet_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeSequenceData)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeSequenceData.subclass:
            return typeSequenceData.subclass(*args_, **kwargs_)
        else:
            return typeSequenceData(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Sequence(self):
        return self.Sequence

    def set_Sequence(self, Sequence):
        self.Sequence = Sequence

    def add_Sequence(self, value):
        self.Sequence.append(value)

    def insert_Sequence_at(self, index, value):
        self.Sequence.insert(index, value)

    def replace_Sequence_at(self, index, value):
        self.Sequence[index] = value

    def get_AuthorSet(self):
        return self.AuthorSet

    def set_AuthorSet(self, AuthorSet):
        self.AuthorSet = AuthorSet

    def hasContent_(self):
        if self.Sequence or self.AuthorSet is not None:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeSequenceData",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeSequenceData")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeSequenceData":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeSequenceData")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeSequenceData", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="typeSequenceData"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeSequenceData",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        for Sequence_ in self.Sequence:
            namespaceprefix_ = self.Sequence_nsprefix_ + ":" if (UseCapturedNS_ and self.Sequence_nsprefix_) else ""
            Sequence_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Sequence", pretty_print=pretty_print
            )
        if self.AuthorSet is not None:
            namespaceprefix_ = self.AuthorSet_nsprefix_ + ":" if (UseCapturedNS_ and self.AuthorSet_nsprefix_) else ""
            self.AuthorSet.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="AuthorSet", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Sequence":
            obj_ = SequenceType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Sequence.append(obj_)
            obj_.original_tagname_ = "Sequence"
        elif nodeName_ == "AuthorSet":
            obj_ = typeAuthorSet.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.AuthorSet = obj_
            obj_.original_tagname_ = "AuthorSet"


# end class typeSequenceData


class typeReleaseStatus(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, Release=None, SetReleaseDate=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.Release = Release
        self.Release_nsprefix_ = None
        self.SetReleaseDate = SetReleaseDate
        self.SetReleaseDate_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeReleaseStatus)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeReleaseStatus.subclass:
            return typeReleaseStatus.subclass(*args_, **kwargs_)
        else:
            return typeReleaseStatus(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Release(self):
        return self.Release

    def set_Release(self, Release):
        self.Release = Release

    def get_SetReleaseDate(self):
        return self.SetReleaseDate

    def set_SetReleaseDate(self, SetReleaseDate):
        self.SetReleaseDate = SetReleaseDate

    def hasContent_(self):
        if self.Release is not None or self.SetReleaseDate is not None:
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="", namespacedef_="", name_="typeReleaseStatus", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeReleaseStatus")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeReleaseStatus":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeReleaseStatus")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile,
                level + 1,
                namespaceprefix_,
                namespacedef_,
                name_="typeReleaseStatus",
                pretty_print=pretty_print,
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="typeReleaseStatus"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="typeReleaseStatus",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Release is not None:
            namespaceprefix_ = self.Release_nsprefix_ + ":" if (UseCapturedNS_ and self.Release_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sRelease>%s</%sRelease>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Release), input_name="Release")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.SetReleaseDate is not None:
            namespaceprefix_ = (
                self.SetReleaseDate_nsprefix_ + ":" if (UseCapturedNS_ and self.SetReleaseDate_nsprefix_) else ""
            )
            self.SetReleaseDate.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="SetReleaseDate", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Release":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Release")
            value_ = self.gds_validate_string(value_, node, "Release")
            self.Release = value_
            self.Release_nsprefix_ = child_.prefix
        elif nodeName_ == "SetReleaseDate":
            obj_ = SetReleaseDateType1.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.SetReleaseDate = obj_
            obj_.original_tagname_ = "SetReleaseDate"


# end class typeReleaseStatus


class Release(GeneratedsSuper):
    """Immediate release"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, Release)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if Release.subclass:
            return Release.subclass(*args_, **kwargs_)
        else:
            return Release(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def hasContent_(self):
        if ():
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="Release", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("Release")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "Release":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="Release")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="Release", pretty_print=pretty_print
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="Release"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="Release",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class Release


class typeLocalId(GeneratedsSuper):
    """Local identifier in submission context. Optional submission id. If
    omitted, the current
    submission is assumed."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, submission_id=None, valueOf_=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.submission_id = _cast(None, submission_id)
        self.submission_id_nsprefix_ = None
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeLocalId)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeLocalId.subclass:
            return typeLocalId.subclass(*args_, **kwargs_)
        else:
            return typeLocalId(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_submission_id(self):
        return self.submission_id

    def set_submission_id(self, submission_id):
        self.submission_id = submission_id

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def hasContent_(self):
        if 1 if type(self.valueOf_) in [int, float] else self.valueOf_:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeLocalId", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeLocalId")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeLocalId":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeLocalId")
        if self.hasContent_():
            outfile.write(">")
            outfile.write(self.convert_unicode(self.valueOf_))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeLocalId", pretty_print=pretty_print
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeLocalId"):
        if self.submission_id is not None and "submission_id" not in already_processed:
            already_processed.add("submission_id")
            outfile.write(
                " submission_id=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(quote_attrib(self.submission_id), input_name="submission_id")
                    ),
                )
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeLocalId",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("submission_id", node)
        if value is not None and "submission_id" not in already_processed:
            already_processed.add("submission_id")
            self.submission_id = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class typeLocalId


class typeSPUID(GeneratedsSuper):
    """Unique identifier in submitter context (Submitter Provided Unique
    ID). DEPRICATED and will be removed: Optional submitter id - eg JGI. If
    omitted, the current submitter is assumed. Will be required: Identifier
    of the submitter namespace: This is a controlled vocabulary on NCBI
    side."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, submitter_id=None, spuid_namespace=None, valueOf_=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.submitter_id = _cast(None, submitter_id)
        self.submitter_id_nsprefix_ = None
        self.spuid_namespace = _cast(None, spuid_namespace)
        self.spuid_namespace_nsprefix_ = None
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeSPUID)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeSPUID.subclass:
            return typeSPUID.subclass(*args_, **kwargs_)
        else:
            return typeSPUID(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_submitter_id(self):
        return self.submitter_id

    def set_submitter_id(self, submitter_id):
        self.submitter_id = submitter_id

    def get_spuid_namespace(self):
        return self.spuid_namespace

    def set_spuid_namespace(self, spuid_namespace):
        self.spuid_namespace = spuid_namespace

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def hasContent_(self):
        if 1 if type(self.valueOf_) in [int, float] else self.valueOf_:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeSPUID", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeSPUID")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeSPUID":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeSPUID")
        if self.hasContent_():
            outfile.write(">")
            outfile.write(self.convert_unicode(self.valueOf_))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeSPUID", pretty_print=pretty_print
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeSPUID"):
        if self.submitter_id is not None and "submitter_id" not in already_processed:
            already_processed.add("submitter_id")
            outfile.write(
                " submitter_id=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.submitter_id), input_name="submitter_id")),)
            )
        if self.spuid_namespace is not None and "spuid_namespace" not in already_processed:
            already_processed.add("spuid_namespace")
            outfile.write(
                " spuid_namespace=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(quote_attrib(self.spuid_namespace), input_name="spuid_namespace")
                    ),
                )
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeSPUID",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("submitter_id", node)
        if value is not None and "submitter_id" not in already_processed:
            already_processed.add("submitter_id")
            self.submitter_id = value
        value = find_attr_value_("spuid_namespace", node)
        if value is not None and "spuid_namespace" not in already_processed:
            already_processed.add("spuid_namespace")
            self.spuid_namespace = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class typeSPUID


class typePrimaryId(GeneratedsSuper):
    """Unique accession in NCBI Archive. Accession is assigned only after
    object is successfully loaded into NCBI archive. Optioanl integer id is for
    internal
    NCBI use only. Optional identifier of the archive. Can be ommitted if
    defined by context. Host archive integer id. May be assigned only by
    NCBI."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, db=None, id=None, valueOf_=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.db = _cast(None, db)
        self.db_nsprefix_ = None
        self.id = _cast(int, id)
        self.id_nsprefix_ = None
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typePrimaryId)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typePrimaryId.subclass:
            return typePrimaryId.subclass(*args_, **kwargs_)
        else:
            return typePrimaryId(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_db(self):
        return self.db

    def set_db(self, db):
        self.db = db

    def get_id(self):
        return self.id

    def set_id(self, id):
        self.id = id

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def hasContent_(self):
        if 1 if type(self.valueOf_) in [int, float] else self.valueOf_:
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typePrimaryId", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typePrimaryId")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typePrimaryId":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typePrimaryId")
        if self.hasContent_():
            outfile.write(">")
            outfile.write(self.convert_unicode(self.valueOf_))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typePrimaryId", pretty_print=pretty_print
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typePrimaryId"):
        if self.db is not None and "db" not in already_processed:
            already_processed.add("db")
            outfile.write(" db=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.db), input_name="db")),))
        if self.id is not None and "id" not in already_processed:
            already_processed.add("id")
            outfile.write(' id="%s"' % self.gds_format_integer(self.id, input_name="id"))

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typePrimaryId",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("db", node)
        if value is not None and "db" not in already_processed:
            already_processed.add("db")
            self.db = value
        value = find_attr_value_("id", node)
        if value is not None and "id" not in already_processed:
            already_processed.add("id")
            self.id = self.gds_parse_integer(value, node, "id")

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class typePrimaryId


class typeIdentifier(GeneratedsSuper):
    """Identifier placed on submitted object by submitter. This is used
    to tie the submitted object with assigned accession."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, PrimaryId=None, SPUID=None, LocalId=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.PrimaryId = PrimaryId
        self.PrimaryId_nsprefix_ = None
        self.SPUID = SPUID
        self.SPUID_nsprefix_ = None
        self.LocalId = LocalId
        self.LocalId_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeIdentifier)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeIdentifier.subclass:
            return typeIdentifier.subclass(*args_, **kwargs_)
        else:
            return typeIdentifier(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_PrimaryId(self):
        return self.PrimaryId

    def set_PrimaryId(self, PrimaryId):
        self.PrimaryId = PrimaryId

    def get_SPUID(self):
        return self.SPUID

    def set_SPUID(self, SPUID):
        self.SPUID = SPUID

    def get_LocalId(self):
        return self.LocalId

    def set_LocalId(self, LocalId):
        self.LocalId = LocalId

    def hasContent_(self):
        if self.PrimaryId is not None or self.SPUID is not None or self.LocalId is not None:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeIdentifier",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeIdentifier")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeIdentifier":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeIdentifier")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeIdentifier", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeIdentifier"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeIdentifier",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.PrimaryId is not None:
            namespaceprefix_ = self.PrimaryId_nsprefix_ + ":" if (UseCapturedNS_ and self.PrimaryId_nsprefix_) else ""
            self.PrimaryId.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="PrimaryId", pretty_print=pretty_print
            )
        if self.SPUID is not None:
            namespaceprefix_ = self.SPUID_nsprefix_ + ":" if (UseCapturedNS_ and self.SPUID_nsprefix_) else ""
            self.SPUID.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="SPUID", pretty_print=pretty_print
            )
        if self.LocalId is not None:
            namespaceprefix_ = self.LocalId_nsprefix_ + ":" if (UseCapturedNS_ and self.LocalId_nsprefix_) else ""
            self.LocalId.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="LocalId", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "PrimaryId":
            obj_ = typePrimaryId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.PrimaryId = obj_
            obj_.original_tagname_ = "PrimaryId"
        elif nodeName_ == "SPUID":
            obj_ = typeSPUID.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.SPUID = obj_
            obj_.original_tagname_ = "SPUID"
        elif nodeName_ == "LocalId":
            obj_ = typeLocalId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.LocalId = obj_
            obj_.original_tagname_ = "LocalId"


# end class typeIdentifier


class typeRefId(GeneratedsSuper):
    """Reference to a record inside NCBI database or in Submission Portal."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, LocalId=None, SPUID=None, PrimaryId=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.LocalId = LocalId
        self.LocalId_nsprefix_ = None
        self.SPUID = SPUID
        self.SPUID_nsprefix_ = None
        self.PrimaryId = PrimaryId
        self.PrimaryId_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeRefId)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeRefId.subclass:
            return typeRefId.subclass(*args_, **kwargs_)
        else:
            return typeRefId(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_LocalId(self):
        return self.LocalId

    def set_LocalId(self, LocalId):
        self.LocalId = LocalId

    def get_SPUID(self):
        return self.SPUID

    def set_SPUID(self, SPUID):
        self.SPUID = SPUID

    def get_PrimaryId(self):
        return self.PrimaryId

    def set_PrimaryId(self, PrimaryId):
        self.PrimaryId = PrimaryId

    def hasContent_(self):
        if self.LocalId is not None or self.SPUID is not None or self.PrimaryId is not None:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeRefId",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeRefId")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeRefId":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeRefId")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeRefId", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeRefId"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeRefId",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.LocalId is not None:
            namespaceprefix_ = self.LocalId_nsprefix_ + ":" if (UseCapturedNS_ and self.LocalId_nsprefix_) else ""
            self.LocalId.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="LocalId", pretty_print=pretty_print
            )
        if self.SPUID is not None:
            namespaceprefix_ = self.SPUID_nsprefix_ + ":" if (UseCapturedNS_ and self.SPUID_nsprefix_) else ""
            self.SPUID.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="SPUID", pretty_print=pretty_print
            )
        if self.PrimaryId is not None:
            namespaceprefix_ = self.PrimaryId_nsprefix_ + ":" if (UseCapturedNS_ and self.PrimaryId_nsprefix_) else ""
            self.PrimaryId.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="PrimaryId", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "LocalId":
            obj_ = typeLocalId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.LocalId = obj_
            obj_.original_tagname_ = "LocalId"
        elif nodeName_ == "SPUID":
            obj_ = typeSPUID.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.SPUID = obj_
            obj_.original_tagname_ = "SPUID"
        elif nodeName_ == "PrimaryId":
            obj_ = typePrimaryId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.PrimaryId = obj_
            obj_.original_tagname_ = "PrimaryId"


# end class typeRefId


class typeLink(GeneratedsSuper):
    """The typeLink represents links between archived objects. If link is
    submitted together with archive object(s) user or local ids can be used.
    Those will
    be resolved once archive objects are accessioned."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, From=None, To=None, Attributes=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.From = From
        self.From_nsprefix_ = None
        self.To = To
        self.To_nsprefix_ = None
        self.Attributes = Attributes
        self.Attributes_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeLink)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeLink.subclass:
            return typeLink.subclass(*args_, **kwargs_)
        else:
            return typeLink(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_From(self):
        return self.From

    def set_From(self, From):
        self.From = From

    def get_To(self):
        return self.To

    def set_To(self, To):
        self.To = To

    def get_Attributes(self):
        return self.Attributes

    def set_Attributes(self, Attributes):
        self.Attributes = Attributes

    def hasContent_(self):
        if self.From is not None or self.To is not None or self.Attributes is not None:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeLink",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeLink")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeLink":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeLink")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeLink", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeLink"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeLink",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.From is not None:
            namespaceprefix_ = self.From_nsprefix_ + ":" if (UseCapturedNS_ and self.From_nsprefix_) else ""
            self.From.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="From", pretty_print=pretty_print
            )
        if self.To is not None:
            namespaceprefix_ = self.To_nsprefix_ + ":" if (UseCapturedNS_ and self.To_nsprefix_) else ""
            self.To.export(outfile, level, namespaceprefix_, namespacedef_="", name_="To", pretty_print=pretty_print)
        if self.Attributes is not None:
            namespaceprefix_ = self.Attributes_nsprefix_ + ":" if (UseCapturedNS_ and self.Attributes_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sAttributes>%s</%sAttributes>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Attributes), input_name="Attributes")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "From":
            obj_ = typeRefId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.From = obj_
            obj_.original_tagname_ = "From"
        elif nodeName_ == "To":
            obj_ = typeRefId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.To = obj_
            obj_.original_tagname_ = "To"
        elif nodeName_ == "Attributes":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Attributes")
            value_ = self.gds_validate_string(value_, node, "Attributes")
            self.Attributes = value_
            self.Attributes_nsprefix_ = child_.prefix


# end class typeLink


class typeExternalLink(GeneratedsSuper):
    """Text which is shown on the obect presentation page for this link."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self, label=None, category=None, URL=None, ExternalID=None, EntrezQuery=None, gds_collector_=None, **kwargs_
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.label = _cast(None, label)
        self.label_nsprefix_ = None
        self.category = _cast(None, category)
        self.category_nsprefix_ = None
        self.URL = URL
        self.URL_nsprefix_ = None
        self.ExternalID = ExternalID
        self.ExternalID_nsprefix_ = None
        self.EntrezQuery = EntrezQuery
        self.EntrezQuery_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeExternalLink)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeExternalLink.subclass:
            return typeExternalLink.subclass(*args_, **kwargs_)
        else:
            return typeExternalLink(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_URL(self):
        return self.URL

    def set_URL(self, URL):
        self.URL = URL

    def get_ExternalID(self):
        return self.ExternalID

    def set_ExternalID(self, ExternalID):
        self.ExternalID = ExternalID

    def get_EntrezQuery(self):
        return self.EntrezQuery

    def set_EntrezQuery(self, EntrezQuery):
        self.EntrezQuery = EntrezQuery

    def get_label(self):
        return self.label

    def set_label(self, label):
        self.label = label

    def get_category(self):
        return self.category

    def set_category(self, category):
        self.category = category

    def hasContent_(self):
        if self.URL is not None or self.ExternalID is not None or self.EntrezQuery is not None:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeExternalLink",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeExternalLink")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeExternalLink":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeExternalLink")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeExternalLink", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeExternalLink"):
        if self.label is not None and "label" not in already_processed:
            already_processed.add("label")
            outfile.write(
                " label=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.label), input_name="label")),)
            )
        if self.category is not None and "category" not in already_processed:
            already_processed.add("category")
            outfile.write(
                " category=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.category), input_name="category")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeExternalLink",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.URL is not None:
            namespaceprefix_ = self.URL_nsprefix_ + ":" if (UseCapturedNS_ and self.URL_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sURL>%s</%sURL>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.URL), input_name="URL")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.ExternalID is not None:
            namespaceprefix_ = self.ExternalID_nsprefix_ + ":" if (UseCapturedNS_ and self.ExternalID_nsprefix_) else ""
            self.ExternalID.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="ExternalID", pretty_print=pretty_print
            )
        if self.EntrezQuery is not None:
            namespaceprefix_ = (
                self.EntrezQuery_nsprefix_ + ":" if (UseCapturedNS_ and self.EntrezQuery_nsprefix_) else ""
            )
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sEntrezQuery>%s</%sEntrezQuery>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.EntrezQuery), input_name="EntrezQuery")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("label", node)
        if value is not None and "label" not in already_processed:
            already_processed.add("label")
            self.label = value
        value = find_attr_value_("category", node)
        if value is not None and "category" not in already_processed:
            already_processed.add("category")
            self.category = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "URL":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "URL")
            value_ = self.gds_validate_string(value_, node, "URL")
            self.URL = value_
            self.URL_nsprefix_ = child_.prefix
        elif nodeName_ == "ExternalID":
            obj_ = typePrimaryId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.ExternalID = obj_
            obj_.original_tagname_ = "ExternalID"
        elif nodeName_ == "EntrezQuery":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "EntrezQuery")
            value_ = self.gds_validate_string(value_, node, "EntrezQuery")
            self.EntrezQuery = value_
            self.EntrezQuery_nsprefix_ = child_.prefix


# end class typeExternalLink


class typeAuthorSet(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, Author=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if Author is None:
            self.Author = []
        else:
            self.Author = Author
        self.Author_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeAuthorSet)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeAuthorSet.subclass:
            return typeAuthorSet.subclass(*args_, **kwargs_)
        else:
            return typeAuthorSet(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Author(self):
        return self.Author

    def set_Author(self, Author):
        self.Author = Author

    def add_Author(self, value):
        self.Author.append(value)

    def insert_Author_at(self, index, value):
        self.Author.insert(index, value)

    def replace_Author_at(self, index, value):
        self.Author[index] = value

    def hasContent_(self):
        if self.Author:
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeAuthorSet", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeAuthorSet")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeAuthorSet":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeAuthorSet")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeAuthorSet", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeAuthorSet"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeAuthorSet",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        for Author_ in self.Author:
            namespaceprefix_ = self.Author_nsprefix_ + ":" if (UseCapturedNS_ and self.Author_nsprefix_) else ""
            Author_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Author", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Author":
            obj_ = AuthorType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Author.append(obj_)
            obj_.original_tagname_ = "Author"


# end class typeAuthorSet


class typePublication(GeneratedsSuper):
    """Unique publication identifier in the specified database that is
    specific to the project. Publication date."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        id=None,
        date=None,
        status=None,
        AuthorSet=None,
        Reference=None,
        StructuredCitation=None,
        DbType=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.id = _cast(None, id)
        self.id_nsprefix_ = None
        if isinstance(date, BaseStrType_):
            initvalue_ = datetime_.datetime.strptime(date, "%Y-%m-%dT%H:%M:%S")
        else:
            initvalue_ = date
        self.date = initvalue_
        self.status = _cast(None, status)
        self.status_nsprefix_ = None
        self.AuthorSet = AuthorSet
        self.AuthorSet_nsprefix_ = None
        self.Reference = Reference
        self.Reference_nsprefix_ = None
        self.StructuredCitation = StructuredCitation
        self.StructuredCitation_nsprefix_ = None
        self.DbType = DbType
        if self.DbType is not None:
            self.validate_DbTypeType(self.DbType)
        self.DbType_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typePublication)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typePublication.subclass:
            return typePublication.subclass(*args_, **kwargs_)
        else:
            return typePublication(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_AuthorSet(self):
        return self.AuthorSet

    def set_AuthorSet(self, AuthorSet):
        self.AuthorSet = AuthorSet

    def get_Reference(self):
        return self.Reference

    def set_Reference(self, Reference):
        self.Reference = Reference

    def get_StructuredCitation(self):
        return self.StructuredCitation

    def set_StructuredCitation(self, StructuredCitation):
        self.StructuredCitation = StructuredCitation

    def get_DbType(self):
        return self.DbType

    def set_DbType(self, DbType):
        self.DbType = DbType

    def get_id(self):
        return self.id

    def set_id(self, id):
        self.id = id

    def get_date(self):
        return self.date

    def set_date(self, date):
        self.date = date

    def get_status(self):
        return self.status

    def set_status(self, status):
        self.status = status

    def validate_DbTypeType(self, value):
        result = True
        # Validate type DbTypeType, a restriction on xs:token.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, str):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (str)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            value = value
            enumerations = ["ePMC", "ePubmed", "eDOI", "eNotAvailable"]
            if value not in enumerations:
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s does not match xsd enumeration restriction on DbTypeType'
                    % {"value": encode_str_2_3(value), "lineno": lineno}
                )
                result = False
        return result

    def validate_statusType(self, value):
        # Validate type statusType, a restriction on xs:token.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, str):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (str)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            value = value
            enumerations = ["ePublished", "eInPress", "eUnpublished"]
            if value not in enumerations:
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s does not match xsd enumeration restriction on statusType'
                    % {"value": encode_str_2_3(value), "lineno": lineno}
                )
                result = False

    def hasContent_(self):
        if (
            self.AuthorSet is not None
            or self.Reference is not None
            or self.StructuredCitation is not None
            or self.DbType is not None
        ):
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typePublication",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typePublication")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typePublication":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typePublication")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typePublication", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typePublication"):
        if self.id is not None and "id" not in already_processed:
            already_processed.add("id")
            outfile.write(" id=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.id), input_name="id")),))
        if self.date is not None and "date" not in already_processed:
            already_processed.add("date")
            outfile.write(' date="%s"' % self.gds_format_datetime(self.date, input_name="date"))
        if self.status is not None and "status" not in already_processed:
            already_processed.add("status")
            outfile.write(
                " status=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.status), input_name="status")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typePublication",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.AuthorSet is not None:
            namespaceprefix_ = self.AuthorSet_nsprefix_ + ":" if (UseCapturedNS_ and self.AuthorSet_nsprefix_) else ""
            self.AuthorSet.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="AuthorSet", pretty_print=pretty_print
            )
        if self.Reference is not None:
            namespaceprefix_ = self.Reference_nsprefix_ + ":" if (UseCapturedNS_ and self.Reference_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sReference>%s</%sReference>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Reference), input_name="Reference")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.StructuredCitation is not None:
            namespaceprefix_ = (
                self.StructuredCitation_nsprefix_ + ":"
                if (UseCapturedNS_ and self.StructuredCitation_nsprefix_)
                else ""
            )
            self.StructuredCitation.export(
                outfile,
                level,
                namespaceprefix_,
                namespacedef_="",
                name_="StructuredCitation",
                pretty_print=pretty_print,
            )
        if self.DbType is not None:
            namespaceprefix_ = self.DbType_nsprefix_ + ":" if (UseCapturedNS_ and self.DbType_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sDbType>%s</%sDbType>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.DbType), input_name="DbType")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("id", node)
        if value is not None and "id" not in already_processed:
            already_processed.add("id")
            self.id = value
        value = find_attr_value_("date", node)
        if value is not None and "date" not in already_processed:
            already_processed.add("date")
            try:
                self.date = self.gds_parse_datetime(value)
            except ValueError as exp:
                raise ValueError("Bad date-time attribute (date): %s" % exp)
        value = find_attr_value_("status", node)
        if value is not None and "status" not in already_processed:
            already_processed.add("status")
            self.status = value
            self.status = " ".join(self.status.split())
            self.validate_statusType(self.status)  # validate type statusType

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "AuthorSet":
            obj_ = typeAuthorSet.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.AuthorSet = obj_
            obj_.original_tagname_ = "AuthorSet"
        elif nodeName_ == "Reference":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Reference")
            value_ = self.gds_validate_string(value_, node, "Reference")
            self.Reference = value_
            self.Reference_nsprefix_ = child_.prefix
        elif nodeName_ == "StructuredCitation":
            obj_ = StructuredCitationType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.StructuredCitation = obj_
            obj_.original_tagname_ = "StructuredCitation"
        elif nodeName_ == "DbType":
            value_ = child_.text
            if value_:
                value_ = re_.sub(String_cleanup_pat_, " ", value_).strip()
            else:
                value_ = ""
            value_ = self.gds_parse_string(value_, node, "DbType")
            value_ = self.gds_validate_string(value_, node, "DbType")
            self.DbType = value_
            self.DbType_nsprefix_ = child_.prefix
            # validate type DbTypeType
            self.validate_DbTypeType(self.DbType)


# end class typePublication


class typeDescriptor(GeneratedsSuper):
    """Common description of an object: title, publication, etc..."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, Title=None, Description=None, ExternalLink=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.Title = Title
        self.Title_nsprefix_ = None
        self.Description = Description
        self.Description_nsprefix_ = None
        if ExternalLink is None:
            self.ExternalLink = []
        else:
            self.ExternalLink = ExternalLink
        self.ExternalLink_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeDescriptor)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeDescriptor.subclass:
            return typeDescriptor.subclass(*args_, **kwargs_)
        else:
            return typeDescriptor(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Title(self):
        return self.Title

    def set_Title(self, Title):
        self.Title = Title

    def get_Description(self):
        return self.Description

    def set_Description(self, Description):
        self.Description = Description

    def get_ExternalLink(self):
        return self.ExternalLink

    def set_ExternalLink(self, ExternalLink):
        self.ExternalLink = ExternalLink

    def add_ExternalLink(self, value):
        self.ExternalLink.append(value)

    def insert_ExternalLink_at(self, index, value):
        self.ExternalLink.insert(index, value)

    def replace_ExternalLink_at(self, index, value):
        self.ExternalLink[index] = value

    def hasContent_(self):
        if self.Title is not None or self.Description is not None or self.ExternalLink:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeDescriptor",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeDescriptor")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeDescriptor":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeDescriptor")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeDescriptor", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeDescriptor"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeDescriptor",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Title is not None:
            namespaceprefix_ = self.Title_nsprefix_ + ":" if (UseCapturedNS_ and self.Title_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sTitle>%s</%sTitle>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Title), input_name="Title")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Description is not None:
            namespaceprefix_ = (
                self.Description_nsprefix_ + ":" if (UseCapturedNS_ and self.Description_nsprefix_) else ""
            )
            self.Description.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Description", pretty_print=pretty_print
            )
        for ExternalLink_ in self.ExternalLink:
            namespaceprefix_ = (
                self.ExternalLink_nsprefix_ + ":" if (UseCapturedNS_ and self.ExternalLink_nsprefix_) else ""
            )
            ExternalLink_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="ExternalLink", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Title":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Title")
            value_ = self.gds_validate_string(value_, node, "Title")
            self.Title = value_
            self.Title_nsprefix_ = child_.prefix
        elif nodeName_ == "Description":
            obj_ = typeBlock.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Description = obj_
            obj_.original_tagname_ = "Description"
        elif nodeName_ == "ExternalLink":
            obj_ = typeExternalLink.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.ExternalLink.append(obj_)
            obj_.original_tagname_ = "ExternalLink"


# end class typeDescriptor


class typeAddress(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        postal_code=None,
        Department=None,
        Institution=None,
        Street=None,
        City=None,
        Sub=None,
        Country=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.postal_code = _cast(None, postal_code)
        self.postal_code_nsprefix_ = None
        self.Department = Department
        self.Department_nsprefix_ = None
        self.Institution = Institution
        self.Institution_nsprefix_ = None
        self.Street = Street
        self.Street_nsprefix_ = None
        self.City = City
        self.City_nsprefix_ = None
        self.Sub = Sub
        self.Sub_nsprefix_ = None
        self.Country = Country
        self.Country_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeAddress)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeAddress.subclass:
            return typeAddress.subclass(*args_, **kwargs_)
        else:
            return typeAddress(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Department(self):
        return self.Department

    def set_Department(self, Department):
        self.Department = Department

    def get_Institution(self):
        return self.Institution

    def set_Institution(self, Institution):
        self.Institution = Institution

    def get_Street(self):
        return self.Street

    def set_Street(self, Street):
        self.Street = Street

    def get_City(self):
        return self.City

    def set_City(self, City):
        self.City = City

    def get_Sub(self):
        return self.Sub

    def set_Sub(self, Sub):
        self.Sub = Sub

    def get_Country(self):
        return self.Country

    def set_Country(self, Country):
        self.Country = Country

    def get_postal_code(self):
        return self.postal_code

    def set_postal_code(self, postal_code):
        self.postal_code = postal_code

    def hasContent_(self):
        if (
            self.Department is not None
            or self.Institution is not None
            or self.Street is not None
            or self.City is not None
            or self.Sub is not None
            or self.Country is not None
        ):
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeAddress", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeAddress")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeAddress":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeAddress")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeAddress", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeAddress"):
        if self.postal_code is not None and "postal_code" not in already_processed:
            already_processed.add("postal_code")
            outfile.write(
                " postal_code=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.postal_code), input_name="postal_code")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeAddress",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Department is not None:
            namespaceprefix_ = self.Department_nsprefix_ + ":" if (UseCapturedNS_ and self.Department_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sDepartment>%s</%sDepartment>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Department), input_name="Department")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Institution is not None:
            namespaceprefix_ = (
                self.Institution_nsprefix_ + ":" if (UseCapturedNS_ and self.Institution_nsprefix_) else ""
            )
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sInstitution>%s</%sInstitution>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Institution), input_name="Institution")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Street is not None:
            namespaceprefix_ = self.Street_nsprefix_ + ":" if (UseCapturedNS_ and self.Street_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sStreet>%s</%sStreet>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Street), input_name="Street")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.City is not None:
            namespaceprefix_ = self.City_nsprefix_ + ":" if (UseCapturedNS_ and self.City_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sCity>%s</%sCity>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.City), input_name="City")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Sub is not None:
            namespaceprefix_ = self.Sub_nsprefix_ + ":" if (UseCapturedNS_ and self.Sub_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sSub>%s</%sSub>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Sub), input_name="Sub")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Country is not None:
            namespaceprefix_ = self.Country_nsprefix_ + ":" if (UseCapturedNS_ and self.Country_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sCountry>%s</%sCountry>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Country), input_name="Country")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("postal_code", node)
        if value is not None and "postal_code" not in already_processed:
            already_processed.add("postal_code")
            self.postal_code = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Department":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Department")
            value_ = self.gds_validate_string(value_, node, "Department")
            self.Department = value_
            self.Department_nsprefix_ = child_.prefix
        elif nodeName_ == "Institution":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Institution")
            value_ = self.gds_validate_string(value_, node, "Institution")
            self.Institution = value_
            self.Institution_nsprefix_ = child_.prefix
        elif nodeName_ == "Street":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Street")
            value_ = self.gds_validate_string(value_, node, "Street")
            self.Street = value_
            self.Street_nsprefix_ = child_.prefix
        elif nodeName_ == "City":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "City")
            value_ = self.gds_validate_string(value_, node, "City")
            self.City = value_
            self.City_nsprefix_ = child_.prefix
        elif nodeName_ == "Sub":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Sub")
            value_ = self.gds_validate_string(value_, node, "Sub")
            self.Sub = value_
            self.Sub_nsprefix_ = child_.prefix
        elif nodeName_ == "Country":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Country")
            value_ = self.gds_validate_string(value_, node, "Country")
            self.Country = value_
            self.Country_nsprefix_ = child_.prefix


# end class typeAddress


class typeName(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, First=None, Last=None, Middle=None, Suffix=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.First = First
        self.First_nsprefix_ = None
        self.Last = Last
        self.Last_nsprefix_ = None
        self.Middle = Middle
        self.Middle_nsprefix_ = None
        self.Suffix = Suffix
        self.Suffix_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeName)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeName.subclass:
            return typeName.subclass(*args_, **kwargs_)
        else:
            return typeName(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_First(self):
        return self.First

    def set_First(self, First):
        self.First = First

    def get_Last(self):
        return self.Last

    def set_Last(self, Last):
        self.Last = Last

    def get_Middle(self):
        return self.Middle

    def set_Middle(self, Middle):
        self.Middle = Middle

    def get_Suffix(self):
        return self.Suffix

    def set_Suffix(self, Suffix):
        self.Suffix = Suffix

    def hasContent_(self):
        if self.First is not None or self.Last is not None or self.Middle is not None or self.Suffix is not None:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeName", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeName")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeName":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeName")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeName", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeName"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeName",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.First is not None:
            namespaceprefix_ = self.First_nsprefix_ + ":" if (UseCapturedNS_ and self.First_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sFirst>%s</%sFirst>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.First), input_name="First")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Last is not None:
            namespaceprefix_ = self.Last_nsprefix_ + ":" if (UseCapturedNS_ and self.Last_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sLast>%s</%sLast>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Last), input_name="Last")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Middle is not None:
            namespaceprefix_ = self.Middle_nsprefix_ + ":" if (UseCapturedNS_ and self.Middle_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sMiddle>%s</%sMiddle>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Middle), input_name="Middle")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Suffix is not None:
            namespaceprefix_ = self.Suffix_nsprefix_ + ":" if (UseCapturedNS_ and self.Suffix_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sSuffix>%s</%sSuffix>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Suffix), input_name="Suffix")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "First":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "First")
            value_ = self.gds_validate_string(value_, node, "First")
            self.First = value_
            self.First_nsprefix_ = child_.prefix
        elif nodeName_ == "Last":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Last")
            value_ = self.gds_validate_string(value_, node, "Last")
            self.Last = value_
            self.Last_nsprefix_ = child_.prefix
        elif nodeName_ == "Middle":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Middle")
            value_ = self.gds_validate_string(value_, node, "Middle")
            self.Middle = value_
            self.Middle_nsprefix_ = child_.prefix
        elif nodeName_ == "Suffix":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Suffix")
            value_ = self.gds_validate_string(value_, node, "Suffix")
            self.Suffix = value_
            self.Suffix_nsprefix_ = child_.prefix


# end class typeName


class typeAuthorName(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, First=None, Last=None, Middle=None, Suffix=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.First = First
        self.First_nsprefix_ = None
        self.Last = Last
        self.Last_nsprefix_ = None
        self.Middle = Middle
        self.Middle_nsprefix_ = None
        self.Suffix = Suffix
        self.Suffix_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeAuthorName)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeAuthorName.subclass:
            return typeAuthorName.subclass(*args_, **kwargs_)
        else:
            return typeAuthorName(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_First(self):
        return self.First

    def set_First(self, First):
        self.First = First

    def get_Last(self):
        return self.Last

    def set_Last(self, Last):
        self.Last = Last

    def get_Middle(self):
        return self.Middle

    def set_Middle(self, Middle):
        self.Middle = Middle

    def get_Suffix(self):
        return self.Suffix

    def set_Suffix(self, Suffix):
        self.Suffix = Suffix

    def hasContent_(self):
        if self.First is not None or self.Last is not None or self.Middle is not None or self.Suffix is not None:
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeAuthorName", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeAuthorName")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeAuthorName":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeAuthorName")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeAuthorName", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeAuthorName"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeAuthorName",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.First is not None:
            namespaceprefix_ = self.First_nsprefix_ + ":" if (UseCapturedNS_ and self.First_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sFirst>%s</%sFirst>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.First), input_name="First")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Last is not None:
            namespaceprefix_ = self.Last_nsprefix_ + ":" if (UseCapturedNS_ and self.Last_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sLast>%s</%sLast>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Last), input_name="Last")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Middle is not None:
            namespaceprefix_ = self.Middle_nsprefix_ + ":" if (UseCapturedNS_ and self.Middle_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sMiddle>%s</%sMiddle>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Middle), input_name="Middle")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Suffix is not None:
            namespaceprefix_ = self.Suffix_nsprefix_ + ":" if (UseCapturedNS_ and self.Suffix_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sSuffix>%s</%sSuffix>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Suffix), input_name="Suffix")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "First":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "First")
            value_ = self.gds_validate_string(value_, node, "First")
            self.First = value_
            self.First_nsprefix_ = child_.prefix
        elif nodeName_ == "Last":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Last")
            value_ = self.gds_validate_string(value_, node, "Last")
            self.Last = value_
            self.Last_nsprefix_ = child_.prefix
        elif nodeName_ == "Middle":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Middle")
            value_ = self.gds_validate_string(value_, node, "Middle")
            self.Middle = value_
            self.Middle_nsprefix_ = child_.prefix
        elif nodeName_ == "Suffix":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Suffix")
            value_ = self.gds_validate_string(value_, node, "Suffix")
            self.Suffix = value_
            self.Suffix_nsprefix_ = child_.prefix


# end class typeAuthorName


class typeContactInfo(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self, email=None, sec_email=None, phone=None, fax=None, Address=None, Name=None, gds_collector_=None, **kwargs_
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.email = _cast(None, email)
        self.email_nsprefix_ = None
        self.sec_email = _cast(None, sec_email)
        self.sec_email_nsprefix_ = None
        self.phone = _cast(None, phone)
        self.phone_nsprefix_ = None
        self.fax = _cast(None, fax)
        self.fax_nsprefix_ = None
        self.Address = Address
        self.Address_nsprefix_ = None
        self.Name = Name
        self.Name_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeContactInfo)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeContactInfo.subclass:
            return typeContactInfo.subclass(*args_, **kwargs_)
        else:
            return typeContactInfo(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Address(self):
        return self.Address

    def set_Address(self, Address):
        self.Address = Address

    def get_Name(self):
        return self.Name

    def set_Name(self, Name):
        self.Name = Name

    def get_email(self):
        return self.email

    def set_email(self, email):
        self.email = email

    def get_sec_email(self):
        return self.sec_email

    def set_sec_email(self, sec_email):
        self.sec_email = sec_email

    def get_phone(self):
        return self.phone

    def set_phone(self, phone):
        self.phone = phone

    def get_fax(self):
        return self.fax

    def set_fax(self, fax):
        self.fax = fax

    def hasContent_(self):
        if self.Address is not None or self.Name is not None:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeContactInfo",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeContactInfo")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeContactInfo":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeContactInfo")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeContactInfo", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeContactInfo"):
        if self.email is not None and "email" not in already_processed:
            already_processed.add("email")
            outfile.write(
                " email=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.email), input_name="email")),)
            )
        if self.sec_email is not None and "sec_email" not in already_processed:
            already_processed.add("sec_email")
            outfile.write(
                " sec_email=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.sec_email), input_name="sec_email")),)
            )
        if self.phone is not None and "phone" not in already_processed:
            already_processed.add("phone")
            outfile.write(
                " phone=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.phone), input_name="phone")),)
            )
        if self.fax is not None and "fax" not in already_processed:
            already_processed.add("fax")
            outfile.write(
                " fax=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.fax), input_name="fax")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeContactInfo",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Address is not None:
            namespaceprefix_ = self.Address_nsprefix_ + ":" if (UseCapturedNS_ and self.Address_nsprefix_) else ""
            self.Address.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Address", pretty_print=pretty_print
            )
        if self.Name is not None:
            namespaceprefix_ = self.Name_nsprefix_ + ":" if (UseCapturedNS_ and self.Name_nsprefix_) else ""
            self.Name.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Name", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("email", node)
        if value is not None and "email" not in already_processed:
            already_processed.add("email")
            self.email = value
        value = find_attr_value_("sec_email", node)
        if value is not None and "sec_email" not in already_processed:
            already_processed.add("sec_email")
            self.sec_email = value
        value = find_attr_value_("phone", node)
        if value is not None and "phone" not in already_processed:
            already_processed.add("phone")
            self.phone = value
        value = find_attr_value_("fax", node)
        if value is not None and "fax" not in already_processed:
            already_processed.add("fax")
            self.fax = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Address":
            obj_ = typeAddress.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Address = obj_
            obj_.original_tagname_ = "Address"
        elif nodeName_ == "Name":
            obj_ = typeName.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Name = obj_
            obj_.original_tagname_ = "Name"


# end class typeContactInfo


class typeOrganism(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        taxonomy_id=None,
        OrganismName=None,
        Label=None,
        Strain=None,
        IsolateName=None,
        Breed=None,
        Cultivar=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.taxonomy_id = _cast(int, taxonomy_id)
        self.taxonomy_id_nsprefix_ = None
        self.OrganismName = OrganismName
        self.OrganismName_nsprefix_ = None
        self.Label = Label
        self.Label_nsprefix_ = None
        self.Strain = Strain
        self.Strain_nsprefix_ = None
        self.IsolateName = IsolateName
        self.IsolateName_nsprefix_ = None
        self.Breed = Breed
        self.Breed_nsprefix_ = None
        self.Cultivar = Cultivar
        self.Cultivar_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeOrganism)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeOrganism.subclass:
            return typeOrganism.subclass(*args_, **kwargs_)
        else:
            return typeOrganism(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_OrganismName(self):
        return self.OrganismName

    def set_OrganismName(self, OrganismName):
        self.OrganismName = OrganismName

    def get_Label(self):
        return self.Label

    def set_Label(self, Label):
        self.Label = Label

    def get_Strain(self):
        return self.Strain

    def set_Strain(self, Strain):
        self.Strain = Strain

    def get_IsolateName(self):
        return self.IsolateName

    def set_IsolateName(self, IsolateName):
        self.IsolateName = IsolateName

    def get_Breed(self):
        return self.Breed

    def set_Breed(self, Breed):
        self.Breed = Breed

    def get_Cultivar(self):
        return self.Cultivar

    def set_Cultivar(self, Cultivar):
        self.Cultivar = Cultivar

    def get_taxonomy_id(self):
        return self.taxonomy_id

    def set_taxonomy_id(self, taxonomy_id):
        self.taxonomy_id = taxonomy_id

    def hasContent_(self):
        if (
            self.OrganismName is not None
            or self.Label is not None
            or self.Strain is not None
            or self.IsolateName is not None
            or self.Breed is not None
            or self.Cultivar is not None
        ):
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeOrganism", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeOrganism")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeOrganism":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeOrganism")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeOrganism", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeOrganism"):
        if self.taxonomy_id is not None and "taxonomy_id" not in already_processed:
            already_processed.add("taxonomy_id")
            outfile.write(' taxonomy_id="%s"' % self.gds_format_integer(self.taxonomy_id, input_name="taxonomy_id"))

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeOrganism",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.OrganismName is not None:
            namespaceprefix_ = (
                self.OrganismName_nsprefix_ + ":" if (UseCapturedNS_ and self.OrganismName_nsprefix_) else ""
            )
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sOrganismName>%s</%sOrganismName>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.OrganismName), input_name="OrganismName")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Label is not None:
            namespaceprefix_ = self.Label_nsprefix_ + ":" if (UseCapturedNS_ and self.Label_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sLabel>%s</%sLabel>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Label), input_name="Label")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Strain is not None:
            namespaceprefix_ = self.Strain_nsprefix_ + ":" if (UseCapturedNS_ and self.Strain_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sStrain>%s</%sStrain>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Strain), input_name="Strain")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.IsolateName is not None:
            namespaceprefix_ = (
                self.IsolateName_nsprefix_ + ":" if (UseCapturedNS_ and self.IsolateName_nsprefix_) else ""
            )
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sIsolateName>%s</%sIsolateName>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.IsolateName), input_name="IsolateName")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Breed is not None:
            namespaceprefix_ = self.Breed_nsprefix_ + ":" if (UseCapturedNS_ and self.Breed_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sBreed>%s</%sBreed>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Breed), input_name="Breed")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Cultivar is not None:
            namespaceprefix_ = self.Cultivar_nsprefix_ + ":" if (UseCapturedNS_ and self.Cultivar_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sCultivar>%s</%sCultivar>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Cultivar), input_name="Cultivar")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("taxonomy_id", node)
        if value is not None and "taxonomy_id" not in already_processed:
            already_processed.add("taxonomy_id")
            self.taxonomy_id = self.gds_parse_integer(value, node, "taxonomy_id")

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "OrganismName":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "OrganismName")
            value_ = self.gds_validate_string(value_, node, "OrganismName")
            self.OrganismName = value_
            self.OrganismName_nsprefix_ = child_.prefix
        elif nodeName_ == "Label":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Label")
            value_ = self.gds_validate_string(value_, node, "Label")
            self.Label = value_
            self.Label_nsprefix_ = child_.prefix
        elif nodeName_ == "Strain":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Strain")
            value_ = self.gds_validate_string(value_, node, "Strain")
            self.Strain = value_
            self.Strain_nsprefix_ = child_.prefix
        elif nodeName_ == "IsolateName":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "IsolateName")
            value_ = self.gds_validate_string(value_, node, "IsolateName")
            self.IsolateName = value_
            self.IsolateName_nsprefix_ = child_.prefix
        elif nodeName_ == "Breed":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Breed")
            value_ = self.gds_validate_string(value_, node, "Breed")
            self.Breed = value_
            self.Breed_nsprefix_ = child_.prefix
        elif nodeName_ == "Cultivar":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Cultivar")
            value_ = self.gds_validate_string(value_, node, "Cultivar")
            self.Cultivar = value_
            self.Cultivar_nsprefix_ = child_.prefix


# end class typeOrganism


class typeBlock(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, p=None, ul=None, ol=None, table=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if p is None:
            self.p = []
        else:
            self.p = p
        self.p_nsprefix_ = None
        if ul is None:
            self.ul = []
        else:
            self.ul = ul
        self.ul_nsprefix_ = None
        if ol is None:
            self.ol = []
        else:
            self.ol = ol
        self.ol_nsprefix_ = None
        if table is None:
            self.table = []
        else:
            self.table = table
        self.table_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeBlock)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeBlock.subclass:
            return typeBlock.subclass(*args_, **kwargs_)
        else:
            return typeBlock(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_p(self):
        return self.p

    def set_p(self, p):
        self.p = p

    def add_p(self, value):
        self.p.append(value)

    def insert_p_at(self, index, value):
        self.p.insert(index, value)

    def replace_p_at(self, index, value):
        self.p[index] = value

    def get_ul(self):
        return self.ul

    def set_ul(self, ul):
        self.ul = ul

    def add_ul(self, value):
        self.ul.append(value)

    def insert_ul_at(self, index, value):
        self.ul.insert(index, value)

    def replace_ul_at(self, index, value):
        self.ul[index] = value

    def get_ol(self):
        return self.ol

    def set_ol(self, ol):
        self.ol = ol

    def add_ol(self, value):
        self.ol.append(value)

    def insert_ol_at(self, index, value):
        self.ol.insert(index, value)

    def replace_ol_at(self, index, value):
        self.ol[index] = value

    def get_table(self):
        return self.table

    def set_table(self, table):
        self.table = table

    def add_table(self, value):
        self.table.append(value)

    def insert_table_at(self, index, value):
        self.table.insert(index, value)

    def replace_table_at(self, index, value):
        self.table[index] = value

    def hasContent_(self):
        if self.p or self.ul or self.ol or self.table:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeBlock",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeBlock")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeBlock":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeBlock")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeBlock", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeBlock"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeBlock",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        for p_ in self.p:
            namespaceprefix_ = self.p_nsprefix_ + ":" if (UseCapturedNS_ and self.p_nsprefix_) else ""
            p_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="p", pretty_print=pretty_print)
        for ul_ in self.ul:
            namespaceprefix_ = self.ul_nsprefix_ + ":" if (UseCapturedNS_ and self.ul_nsprefix_) else ""
            ul_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="ul", pretty_print=pretty_print)
        for ol_ in self.ol:
            namespaceprefix_ = self.ol_nsprefix_ + ":" if (UseCapturedNS_ and self.ol_nsprefix_) else ""
            ol_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="ol", pretty_print=pretty_print)
        for table_ in self.table:
            namespaceprefix_ = self.table_nsprefix_ + ":" if (UseCapturedNS_ and self.table_nsprefix_) else ""
            table_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="table", pretty_print=pretty_print)

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "p":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            obj_ = class_obj_.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.p.append(obj_)
            obj_.original_tagname_ = "p"
        elif nodeName_ == "ul":
            obj_ = typeL.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.ul.append(obj_)
            obj_.original_tagname_ = "ul"
        elif nodeName_ == "ol":
            obj_ = typeL.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.ol.append(obj_)
            obj_.original_tagname_ = "ol"
        elif nodeName_ == "table":
            obj_ = typeTable.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.table.append(obj_)
            obj_.original_tagname_ = "table"


# end class typeBlock


class typeInline(GeneratedsSuper):
    """ "Inline" covers inline or "text-level" elements"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        a=None,
        i=None,
        b=None,
        sub=None,
        sup=None,
        valueOf_=None,
        mixedclass_=None,
        content_=None,
        extensiontype_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if a is None:
            self.a = []
        else:
            self.a = a
        self.a_nsprefix_ = None
        if i is None:
            self.i = []
        else:
            self.i = i
        self.i_nsprefix_ = None
        if b is None:
            self.b = []
        else:
            self.b = b
        self.b_nsprefix_ = None
        if sub is None:
            self.sub = []
        else:
            self.sub = sub
        self.sub_nsprefix_ = None
        if sup is None:
            self.sup = []
        else:
            self.sup = sup
        self.sup_nsprefix_ = None
        self.valueOf_ = valueOf_
        self.extensiontype_ = extensiontype_
        if mixedclass_ is None:
            self.mixedclass_ = MixedContainer
        else:
            self.mixedclass_ = mixedclass_
        if content_ is None:
            self.content_ = []
        else:
            self.content_ = content_
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeInline)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeInline.subclass:
            return typeInline.subclass(*args_, **kwargs_)
        else:
            return typeInline(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_a(self):
        return self.a

    def set_a(self, a):
        self.a = a

    def add_a(self, value):
        self.a.append(value)

    def insert_a_at(self, index, value):
        self.a.insert(index, value)

    def replace_a_at(self, index, value):
        self.a[index] = value

    def get_i(self):
        return self.i

    def set_i(self, i):
        self.i = i

    def add_i(self, value):
        self.i.append(value)

    def insert_i_at(self, index, value):
        self.i.insert(index, value)

    def replace_i_at(self, index, value):
        self.i[index] = value

    def get_b(self):
        return self.b

    def set_b(self, b):
        self.b = b

    def add_b(self, value):
        self.b.append(value)

    def insert_b_at(self, index, value):
        self.b.insert(index, value)

    def replace_b_at(self, index, value):
        self.b[index] = value

    def get_sub(self):
        return self.sub

    def set_sub(self, sub):
        self.sub = sub

    def add_sub(self, value):
        self.sub.append(value)

    def insert_sub_at(self, index, value):
        self.sub.insert(index, value)

    def replace_sub_at(self, index, value):
        self.sub[index] = value

    def get_sup(self):
        return self.sup

    def set_sup(self, sup):
        self.sup = sup

    def add_sup(self, value):
        self.sup.append(value)

    def insert_sup_at(self, index, value):
        self.sup.insert(index, value)

    def replace_sup_at(self, index, value):
        self.sup[index] = value

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def get_extensiontype_(self):
        return self.extensiontype_

    def set_extensiontype_(self, extensiontype_):
        self.extensiontype_ = extensiontype_

    def hasContent_(self):
        if (
            self.a
            or self.i
            or self.b
            or self.sub
            or self.sup
            or (1 if type(self.valueOf_) in [int, float] else self.valueOf_)
            or self.content_
        ):
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeInline",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeInline")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeInline":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeInline")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeInline", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeInline"):
        if self.extensiontype_ is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            outfile.write(' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"')
            if ":" not in self.extensiontype_:
                imported_ns_type_prefix_ = GenerateDSNamespaceTypePrefixes_.get(self.extensiontype_, "")
                outfile.write(' xsi:type="%s%s"' % (imported_ns_type_prefix_, self.extensiontype_))
            else:
                outfile.write(' xsi:type="%s"' % self.extensiontype_)
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeInline",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if not fromsubclass_:
            for item_ in self.content_:
                item_.export(outfile, level, item_.name, namespaceprefix_, pretty_print=pretty_print)
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        for a_ in self.a:
            namespaceprefix_ = self.a_nsprefix_ + ":" if (UseCapturedNS_ and self.a_nsprefix_) else ""
            a_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="a", pretty_print=pretty_print)
        for i_ in self.i:
            namespaceprefix_ = self.i_nsprefix_ + ":" if (UseCapturedNS_ and self.i_nsprefix_) else ""
            i_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="i", pretty_print=pretty_print)
        for b_ in self.b:
            namespaceprefix_ = self.b_nsprefix_ + ":" if (UseCapturedNS_ and self.b_nsprefix_) else ""
            b_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="b", pretty_print=pretty_print)
        for sub_ in self.sub:
            namespaceprefix_ = self.sub_nsprefix_ + ":" if (UseCapturedNS_ and self.sub_nsprefix_) else ""
            sub_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="sub", pretty_print=pretty_print)
        for sup_ in self.sup:
            namespaceprefix_ = self.sup_nsprefix_ + ":" if (UseCapturedNS_ and self.sup_nsprefix_) else ""
            sup_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="sup", pretty_print=pretty_print)

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        if node.text is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", node.text)
            self.content_.append(obj_)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("xsi:type", node)
        if value is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            self.extensiontype_ = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "a":
            obj_ = typeA.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "a", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_a"):
                self.add_a(obj_.value)
            elif hasattr(self, "set_a"):
                self.set_a(obj_.value)
        elif nodeName_ == "i":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "i", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_i"):
                self.add_i(obj_.value)
            elif hasattr(self, "set_i"):
                self.set_i(obj_.value)
        elif nodeName_ == "b":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "b", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_b"):
                self.add_b(obj_.value)
            elif hasattr(self, "set_b"):
                self.set_b(obj_.value)
        elif nodeName_ == "sub":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "sub", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_sub"):
                self.add_sub(obj_.value)
            elif hasattr(self, "set_sub"):
                self.set_sub(obj_.value)
        elif nodeName_ == "sup":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "sup", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_sup"):
                self.add_sup(obj_.value)
            elif hasattr(self, "set_sup"):
                self.set_sup(obj_.value)
        if not fromsubclass_ and child_.tail is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", child_.tail)
            self.content_.append(obj_)


# end class typeInline


class typeFlow(GeneratedsSuper):
    """ "Flow" mixes block and inline and is used for list items etc."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        p=None,
        ul=None,
        ol=None,
        table=None,
        a=None,
        i=None,
        b=None,
        sub=None,
        sup=None,
        valueOf_=None,
        mixedclass_=None,
        content_=None,
        extensiontype_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.p = p
        self.p_nsprefix_ = None
        self.ul = ul
        self.ul_nsprefix_ = None
        self.ol = ol
        self.ol_nsprefix_ = None
        self.table = table
        self.table_nsprefix_ = None
        self.a = a
        self.a_nsprefix_ = None
        self.i = i
        self.i_nsprefix_ = None
        self.b = b
        self.b_nsprefix_ = None
        self.sub = sub
        self.sub_nsprefix_ = None
        self.sup = sup
        self.sup_nsprefix_ = None
        self.valueOf_ = valueOf_
        self.extensiontype_ = extensiontype_
        if mixedclass_ is None:
            self.mixedclass_ = MixedContainer
        else:
            self.mixedclass_ = mixedclass_
        if content_ is None:
            self.content_ = []
        else:
            self.content_ = content_
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeFlow)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeFlow.subclass:
            return typeFlow.subclass(*args_, **kwargs_)
        else:
            return typeFlow(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_p(self):
        return self.p

    def set_p(self, p):
        self.p = p

    def get_ul(self):
        return self.ul

    def set_ul(self, ul):
        self.ul = ul

    def get_ol(self):
        return self.ol

    def set_ol(self, ol):
        self.ol = ol

    def get_table(self):
        return self.table

    def set_table(self, table):
        self.table = table

    def get_a(self):
        return self.a

    def set_a(self, a):
        self.a = a

    def get_i(self):
        return self.i

    def set_i(self, i):
        self.i = i

    def get_b(self):
        return self.b

    def set_b(self, b):
        self.b = b

    def get_sub(self):
        return self.sub

    def set_sub(self, sub):
        self.sub = sub

    def get_sup(self):
        return self.sup

    def set_sup(self, sup):
        self.sup = sup

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def get_extensiontype_(self):
        return self.extensiontype_

    def set_extensiontype_(self, extensiontype_):
        self.extensiontype_ = extensiontype_

    def hasContent_(self):
        if (
            self.p is not None
            or self.ul is not None
            or self.ol is not None
            or self.table is not None
            or self.a is not None
            or self.i is not None
            or self.b is not None
            or self.sub is not None
            or self.sup is not None
            or (1 if type(self.valueOf_) in [int, float] else self.valueOf_)
            or self.content_
        ):
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeFlow",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeFlow")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeFlow":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeFlow")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeFlow", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeFlow"):
        if self.extensiontype_ is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            outfile.write(' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"')
            if ":" not in self.extensiontype_:
                imported_ns_type_prefix_ = GenerateDSNamespaceTypePrefixes_.get(self.extensiontype_, "")
                outfile.write(' xsi:type="%s%s"' % (imported_ns_type_prefix_, self.extensiontype_))
            else:
                outfile.write(' xsi:type="%s"' % self.extensiontype_)
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeFlow",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if not fromsubclass_:
            for item_ in self.content_:
                item_.export(outfile, level, item_.name, namespaceprefix_, pretty_print=pretty_print)
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.p is not None:
            namespaceprefix_ = self.p_nsprefix_ + ":" if (UseCapturedNS_ and self.p_nsprefix_) else ""
            self.p.export(outfile, level, namespaceprefix_, namespacedef_="", name_="p", pretty_print=pretty_print)
        if self.ul is not None:
            namespaceprefix_ = self.ul_nsprefix_ + ":" if (UseCapturedNS_ and self.ul_nsprefix_) else ""
            self.ul.export(outfile, level, namespaceprefix_, namespacedef_="", name_="ul", pretty_print=pretty_print)
        if self.ol is not None:
            namespaceprefix_ = self.ol_nsprefix_ + ":" if (UseCapturedNS_ and self.ol_nsprefix_) else ""
            self.ol.export(outfile, level, namespaceprefix_, namespacedef_="", name_="ol", pretty_print=pretty_print)
        if self.table is not None:
            namespaceprefix_ = self.table_nsprefix_ + ":" if (UseCapturedNS_ and self.table_nsprefix_) else ""
            self.table.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="table", pretty_print=pretty_print
            )
        if self.a is not None:
            namespaceprefix_ = self.a_nsprefix_ + ":" if (UseCapturedNS_ and self.a_nsprefix_) else ""
            self.a.export(outfile, level, namespaceprefix_, namespacedef_="", name_="a", pretty_print=pretty_print)
        if self.i is not None:
            namespaceprefix_ = self.i_nsprefix_ + ":" if (UseCapturedNS_ and self.i_nsprefix_) else ""
            self.i.export(outfile, level, namespaceprefix_, namespacedef_="", name_="i", pretty_print=pretty_print)
        if self.b is not None:
            namespaceprefix_ = self.b_nsprefix_ + ":" if (UseCapturedNS_ and self.b_nsprefix_) else ""
            self.b.export(outfile, level, namespaceprefix_, namespacedef_="", name_="b", pretty_print=pretty_print)
        if self.sub is not None:
            namespaceprefix_ = self.sub_nsprefix_ + ":" if (UseCapturedNS_ and self.sub_nsprefix_) else ""
            self.sub.export(outfile, level, namespaceprefix_, namespacedef_="", name_="sub", pretty_print=pretty_print)
        if self.sup is not None:
            namespaceprefix_ = self.sup_nsprefix_ + ":" if (UseCapturedNS_ and self.sup_nsprefix_) else ""
            self.sup.export(outfile, level, namespaceprefix_, namespacedef_="", name_="sup", pretty_print=pretty_print)

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        if node.text is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", node.text)
            self.content_.append(obj_)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("xsi:type", node)
        if value is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            self.extensiontype_ = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "p":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "p", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_p"):
                self.add_p(obj_.value)
            elif hasattr(self, "set_p"):
                self.set_p(obj_.value)
        elif nodeName_ == "ul":
            obj_ = typeL.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "ul", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_ul"):
                self.add_ul(obj_.value)
            elif hasattr(self, "set_ul"):
                self.set_ul(obj_.value)
        elif nodeName_ == "ol":
            obj_ = typeL.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "ol", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_ol"):
                self.add_ol(obj_.value)
            elif hasattr(self, "set_ol"):
                self.set_ol(obj_.value)
        elif nodeName_ == "table":
            obj_ = typeTable.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "table", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_table"):
                self.add_table(obj_.value)
            elif hasattr(self, "set_table"):
                self.set_table(obj_.value)
        elif nodeName_ == "a":
            obj_ = typeA.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "a", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_a"):
                self.add_a(obj_.value)
            elif hasattr(self, "set_a"):
                self.set_a(obj_.value)
        elif nodeName_ == "i":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "i", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_i"):
                self.add_i(obj_.value)
            elif hasattr(self, "set_i"):
                self.set_i(obj_.value)
        elif nodeName_ == "b":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "b", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_b"):
                self.add_b(obj_.value)
            elif hasattr(self, "set_b"):
                self.set_b(obj_.value)
        elif nodeName_ == "sub":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "sub", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_sub"):
                self.add_sub(obj_.value)
            elif hasattr(self, "set_sub"):
                self.set_sub(obj_.value)
        elif nodeName_ == "sup":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "sup", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_sup"):
                self.add_sup(obj_.value)
            elif hasattr(self, "set_sup"):
                self.set_sup(obj_.value)
        if not fromsubclass_ and child_.tail is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", child_.tail)
            self.content_.append(obj_)


# end class typeFlow


class typeA_content(GeneratedsSuper):
    """a elements use "Inline" excluding a"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        i=None,
        b=None,
        sub=None,
        sup=None,
        valueOf_=None,
        mixedclass_=None,
        content_=None,
        extensiontype_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if i is None:
            self.i = []
        else:
            self.i = i
        self.i_nsprefix_ = None
        if b is None:
            self.b = []
        else:
            self.b = b
        self.b_nsprefix_ = None
        if sub is None:
            self.sub = []
        else:
            self.sub = sub
        self.sub_nsprefix_ = None
        if sup is None:
            self.sup = []
        else:
            self.sup = sup
        self.sup_nsprefix_ = None
        self.valueOf_ = valueOf_
        self.extensiontype_ = extensiontype_
        if mixedclass_ is None:
            self.mixedclass_ = MixedContainer
        else:
            self.mixedclass_ = mixedclass_
        if content_ is None:
            self.content_ = []
        else:
            self.content_ = content_
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeA_content)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeA_content.subclass:
            return typeA_content.subclass(*args_, **kwargs_)
        else:
            return typeA_content(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_i(self):
        return self.i

    def set_i(self, i):
        self.i = i

    def add_i(self, value):
        self.i.append(value)

    def insert_i_at(self, index, value):
        self.i.insert(index, value)

    def replace_i_at(self, index, value):
        self.i[index] = value

    def get_b(self):
        return self.b

    def set_b(self, b):
        self.b = b

    def add_b(self, value):
        self.b.append(value)

    def insert_b_at(self, index, value):
        self.b.insert(index, value)

    def replace_b_at(self, index, value):
        self.b[index] = value

    def get_sub(self):
        return self.sub

    def set_sub(self, sub):
        self.sub = sub

    def add_sub(self, value):
        self.sub.append(value)

    def insert_sub_at(self, index, value):
        self.sub.insert(index, value)

    def replace_sub_at(self, index, value):
        self.sub[index] = value

    def get_sup(self):
        return self.sup

    def set_sup(self, sup):
        self.sup = sup

    def add_sup(self, value):
        self.sup.append(value)

    def insert_sup_at(self, index, value):
        self.sup.insert(index, value)

    def replace_sup_at(self, index, value):
        self.sup[index] = value

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def get_extensiontype_(self):
        return self.extensiontype_

    def set_extensiontype_(self, extensiontype_):
        self.extensiontype_ = extensiontype_

    def hasContent_(self):
        if (
            self.i
            or self.b
            or self.sub
            or self.sup
            or (1 if type(self.valueOf_) in [int, float] else self.valueOf_)
            or self.content_
        ):
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeA.content",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeA.content")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeA.content":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeA.content")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeA.content", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeA.content"):
        if self.extensiontype_ is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            outfile.write(' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"')
            if ":" not in self.extensiontype_:
                imported_ns_type_prefix_ = GenerateDSNamespaceTypePrefixes_.get(self.extensiontype_, "")
                outfile.write(' xsi:type="%s%s"' % (imported_ns_type_prefix_, self.extensiontype_))
            else:
                outfile.write(' xsi:type="%s"' % self.extensiontype_)
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeA.content",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if not fromsubclass_:
            for item_ in self.content_:
                item_.export(outfile, level, item_.name, namespaceprefix_, pretty_print=pretty_print)
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        for i_ in self.i:
            namespaceprefix_ = self.i_nsprefix_ + ":" if (UseCapturedNS_ and self.i_nsprefix_) else ""
            i_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="i", pretty_print=pretty_print)
        for b_ in self.b:
            namespaceprefix_ = self.b_nsprefix_ + ":" if (UseCapturedNS_ and self.b_nsprefix_) else ""
            b_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="b", pretty_print=pretty_print)
        for sub_ in self.sub:
            namespaceprefix_ = self.sub_nsprefix_ + ":" if (UseCapturedNS_ and self.sub_nsprefix_) else ""
            sub_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="sub", pretty_print=pretty_print)
        for sup_ in self.sup:
            namespaceprefix_ = self.sup_nsprefix_ + ":" if (UseCapturedNS_ and self.sup_nsprefix_) else ""
            sup_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="sup", pretty_print=pretty_print)

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        if node.text is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", node.text)
            self.content_.append(obj_)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("xsi:type", node)
        if value is not None and "xsi:type" not in already_processed:
            already_processed.add("xsi:type")
            self.extensiontype_ = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "i":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "i", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_i"):
                self.add_i(obj_.value)
            elif hasattr(self, "set_i"):
                self.set_i(obj_.value)
        elif nodeName_ == "b":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "b", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_b"):
                self.add_b(obj_.value)
            elif hasattr(self, "set_b"):
                self.set_b(obj_.value)
        elif nodeName_ == "sub":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "sub", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_sub"):
                self.add_sub(obj_.value)
            elif hasattr(self, "set_sub"):
                self.set_sub(obj_.value)
        elif nodeName_ == "sup":
            class_obj_ = self.get_class_obj_(child_, typeInline)
            class_obj_ = typeInline.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            obj_ = self.mixedclass_(MixedContainer.CategoryComplex, MixedContainer.TypeNone, "sup", obj_)
            self.content_.append(obj_)
            if hasattr(self, "add_sup"):
                self.add_sup(obj_.value)
            elif hasattr(self, "set_sup"):
                self.set_sup(obj_.value)
        if not fromsubclass_ and child_.tail is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", child_.tail)
            self.content_.append(obj_)


# end class typeA_content


class typeL(GeneratedsSuper):
    """Ordered or Unordered list"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, li=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if li is None:
            self.li = []
        else:
            self.li = li
        self.li_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeL)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeL.subclass:
            return typeL.subclass(*args_, **kwargs_)
        else:
            return typeL(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_li(self):
        return self.li

    def set_li(self, li):
        self.li = li

    def add_li(self, value):
        self.li.append(value)

    def insert_li_at(self, index, value):
        self.li.insert(index, value)

    def replace_li_at(self, index, value):
        self.li[index] = value

    def hasContent_(self):
        if self.li:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeL",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeL")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeL":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeL")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeL", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeL"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeL",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        for li_ in self.li:
            namespaceprefix_ = self.li_nsprefix_ + ":" if (UseCapturedNS_ and self.li_nsprefix_) else ""
            li_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="li", pretty_print=pretty_print)

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "li":
            obj_ = typeLI.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.li.append(obj_)
            obj_.original_tagname_ = "li"


# end class typeL


class typeLI(typeFlow):
    """list item"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = typeFlow

    def __init__(
        self,
        p=None,
        ul=None,
        ol=None,
        table=None,
        a=None,
        i=None,
        b=None,
        sub=None,
        sup=None,
        valueOf_=None,
        mixedclass_=None,
        content_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        super(typeLI, self).__init__(p, ul, ol, table, a, i, b, sub, sup, valueOf_, mixedclass_, content_, **kwargs_)
        self.valueOf_ = valueOf_
        if mixedclass_ is None:
            self.mixedclass_ = MixedContainer
        else:
            self.mixedclass_ = mixedclass_
        if content_ is None:
            self.content_ = []
        else:
            self.content_ = content_
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeLI)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeLI.subclass:
            return typeLI.subclass(*args_, **kwargs_)
        else:
            return typeLI(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def hasContent_(self):
        if (
            (1 if type(self.valueOf_) in [int, float] else self.valueOf_)
            or self.content_
            or super(typeLI, self).hasContent_()
        ):
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeLI", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeLI")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeLI":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeLI")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeLI", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeLI"):
        super(typeLI, self).exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeLI")

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeLI",
        fromsubclass_=False,
        pretty_print=True,
    ):
        super(typeLI, self).exportChildren(
            outfile, level, namespaceprefix_, namespacedef_, name_, True, pretty_print=pretty_print
        )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        if node.text is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", node.text)
            self.content_.append(obj_)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        super(typeLI, self).buildAttributes(node, attrs, already_processed)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if not fromsubclass_ and child_.tail is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", child_.tail)
            self.content_.append(obj_)
        super(typeLI, self).buildChildren(child_, node, nodeName_, True)
        pass


# end class typeLI


class typeA(typeA_content):
    """content is "Inline" except that anchors shouldn't be nested"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = typeA_content

    def __init__(
        self,
        i=None,
        b=None,
        sub=None,
        sup=None,
        href=None,
        type_=None,
        valueOf_=None,
        mixedclass_=None,
        content_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        super(typeA, self).__init__(i, b, sub, sup, valueOf_, mixedclass_, content_, **kwargs_)
        self.href = _cast(None, href)
        self.href_nsprefix_ = None
        self.type_ = _cast(None, type_)
        self.type__nsprefix_ = None
        self.valueOf_ = valueOf_
        if mixedclass_ is None:
            self.mixedclass_ = MixedContainer
        else:
            self.mixedclass_ = mixedclass_
        if content_ is None:
            self.content_ = []
        else:
            self.content_ = content_
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeA)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeA.subclass:
            return typeA.subclass(*args_, **kwargs_)
        else:
            return typeA(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_href(self):
        return self.href

    def set_href(self, href):
        self.href = href

    def get_type(self):
        return self.type_

    def set_type(self, type_):
        self.type_ = type_

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def hasContent_(self):
        if (
            (1 if type(self.valueOf_) in [int, float] else self.valueOf_)
            or self.content_
            or super(typeA, self).hasContent_()
        ):
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeA", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeA")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeA":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeA")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeA", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeA"):
        super(typeA, self).exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeA")
        if self.href is not None and "href" not in already_processed:
            already_processed.add("href")
            outfile.write(
                " href=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.href), input_name="href")),)
            )
        if self.type_ is not None and "type_" not in already_processed:
            already_processed.add("type_")
            outfile.write(
                " type=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.type_), input_name="type")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeA",
        fromsubclass_=False,
        pretty_print=True,
    ):
        super(typeA, self).exportChildren(
            outfile, level, namespaceprefix_, namespacedef_, name_, True, pretty_print=pretty_print
        )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        if node.text is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", node.text)
            self.content_.append(obj_)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("href", node)
        if value is not None and "href" not in already_processed:
            already_processed.add("href")
            self.href = value
        value = find_attr_value_("type", node)
        if value is not None and "type" not in already_processed:
            already_processed.add("type")
            self.type_ = value
        super(typeA, self).buildAttributes(node, attrs, already_processed)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if not fromsubclass_ and child_.tail is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", child_.tail)
            self.content_.append(obj_)
        super(typeA, self).buildChildren(child_, node, nodeName_, True)
        pass


# end class typeA


class typeTable(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, caption=None, tr=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.caption = caption
        self.caption_nsprefix_ = None
        if tr is None:
            self.tr = []
        else:
            self.tr = tr
        self.tr_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeTable)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeTable.subclass:
            return typeTable.subclass(*args_, **kwargs_)
        else:
            return typeTable(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_caption(self):
        return self.caption

    def set_caption(self, caption):
        self.caption = caption

    def get_tr(self):
        return self.tr

    def set_tr(self, tr):
        self.tr = tr

    def add_tr(self, value):
        self.tr.append(value)

    def insert_tr_at(self, index, value):
        self.tr.insert(index, value)

    def replace_tr_at(self, index, value):
        self.tr[index] = value

    def hasContent_(self):
        if self.caption is not None or self.tr:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeTable",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeTable")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeTable":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeTable")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeTable", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeTable"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeTable",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.caption is not None:
            namespaceprefix_ = self.caption_nsprefix_ + ":" if (UseCapturedNS_ and self.caption_nsprefix_) else ""
            self.caption.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="caption", pretty_print=pretty_print
            )
        for tr_ in self.tr:
            namespaceprefix_ = self.tr_nsprefix_ + ":" if (UseCapturedNS_ and self.tr_nsprefix_) else ""
            tr_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="tr", pretty_print=pretty_print)

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "caption":
            obj_ = typeCaption.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.caption = obj_
            obj_.original_tagname_ = "caption"
        elif nodeName_ == "tr":
            obj_ = typeTR.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.tr.append(obj_)
            obj_.original_tagname_ = "tr"


# end class typeTable


class typeCaption(typeInline):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = typeInline

    def __init__(
        self,
        a=None,
        i=None,
        b=None,
        sub=None,
        sup=None,
        valueOf_=None,
        mixedclass_=None,
        content_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        super(typeCaption, self).__init__(a, i, b, sub, sup, valueOf_, mixedclass_, content_, **kwargs_)
        self.valueOf_ = valueOf_
        if mixedclass_ is None:
            self.mixedclass_ = MixedContainer
        else:
            self.mixedclass_ = mixedclass_
        if content_ is None:
            self.content_ = []
        else:
            self.content_ = content_
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeCaption)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeCaption.subclass:
            return typeCaption.subclass(*args_, **kwargs_)
        else:
            return typeCaption(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def hasContent_(self):
        if (
            (1 if type(self.valueOf_) in [int, float] else self.valueOf_)
            or self.content_
            or super(typeCaption, self).hasContent_()
        ):
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeCaption", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeCaption")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeCaption":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeCaption")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeCaption", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeCaption"):
        super(typeCaption, self).exportAttributes(
            outfile, level, already_processed, namespaceprefix_, name_="typeCaption"
        )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeCaption",
        fromsubclass_=False,
        pretty_print=True,
    ):
        super(typeCaption, self).exportChildren(
            outfile, level, namespaceprefix_, namespacedef_, name_, True, pretty_print=pretty_print
        )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        if node.text is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", node.text)
            self.content_.append(obj_)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        super(typeCaption, self).buildAttributes(node, attrs, already_processed)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if not fromsubclass_ and child_.tail is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", child_.tail)
            self.content_.append(obj_)
        super(typeCaption, self).buildChildren(child_, node, nodeName_, True)
        pass


# end class typeCaption


class typeTR(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, th=None, td=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if th is None:
            self.th = []
        else:
            self.th = th
        self.th_nsprefix_ = None
        if td is None:
            self.td = []
        else:
            self.td = td
        self.td_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeTR)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeTR.subclass:
            return typeTR.subclass(*args_, **kwargs_)
        else:
            return typeTR(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_th(self):
        return self.th

    def set_th(self, th):
        self.th = th

    def add_th(self, value):
        self.th.append(value)

    def insert_th_at(self, index, value):
        self.th.insert(index, value)

    def replace_th_at(self, index, value):
        self.th[index] = value

    def get_td(self):
        return self.td

    def set_td(self, td):
        self.td = td

    def add_td(self, value):
        self.td.append(value)

    def insert_td_at(self, index, value):
        self.td.insert(index, value)

    def replace_td_at(self, index, value):
        self.td[index] = value

    def hasContent_(self):
        if self.th or self.td:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeTR",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeTR")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeTR":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeTR")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeTR", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeTR"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="typeTR",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        for th_ in self.th:
            namespaceprefix_ = self.th_nsprefix_ + ":" if (UseCapturedNS_ and self.th_nsprefix_) else ""
            th_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="th", pretty_print=pretty_print)
        for td_ in self.td:
            namespaceprefix_ = self.td_nsprefix_ + ":" if (UseCapturedNS_ and self.td_nsprefix_) else ""
            td_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="td", pretty_print=pretty_print)

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "th":
            obj_ = typeTH.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.th.append(obj_)
            obj_.original_tagname_ = "th"
        elif nodeName_ == "td":
            obj_ = typeTD.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.td.append(obj_)
            obj_.original_tagname_ = "td"


# end class typeTR


class typeTH(typeFlow):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = typeFlow

    def __init__(
        self,
        p=None,
        ul=None,
        ol=None,
        table=None,
        a=None,
        i=None,
        b=None,
        sub=None,
        sup=None,
        rowspan="1",
        colspan="1",
        valueOf_=None,
        mixedclass_=None,
        content_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        super(typeTH, self).__init__(p, ul, ol, table, a, i, b, sub, sup, valueOf_, mixedclass_, content_, **kwargs_)
        self.rowspan = _cast(None, rowspan)
        self.rowspan_nsprefix_ = None
        self.colspan = _cast(None, colspan)
        self.colspan_nsprefix_ = None
        self.valueOf_ = valueOf_
        if mixedclass_ is None:
            self.mixedclass_ = MixedContainer
        else:
            self.mixedclass_ = mixedclass_
        if content_ is None:
            self.content_ = []
        else:
            self.content_ = content_
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeTH)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeTH.subclass:
            return typeTH.subclass(*args_, **kwargs_)
        else:
            return typeTH(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_rowspan(self):
        return self.rowspan

    def set_rowspan(self, rowspan):
        self.rowspan = rowspan

    def get_colspan(self):
        return self.colspan

    def set_colspan(self, colspan):
        self.colspan = colspan

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def validate_typeNumber(self, value):
        # Validate type com:typeNumber, a restriction on xs:nonNegativeInteger.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, int):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (int)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            if not self.gds_validate_simple_patterns(self.validate_typeNumber_patterns_, value):
                self.gds_collector_.add_message(
                    'Value "%s" does not match xsd pattern restrictions: %s'
                    % (
                        encode_str_2_3(value),
                        self.validate_typeNumber_patterns_,
                    )
                )

    validate_typeNumber_patterns_ = [["^([0-9]+)$"]]

    def hasContent_(self):
        if (
            (1 if type(self.valueOf_) in [int, float] else self.valueOf_)
            or self.content_
            or super(typeTH, self).hasContent_()
        ):
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeTH", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeTH")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeTH":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeTH")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeTH", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeTH"):
        super(typeTH, self).exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeTH")
        if self.rowspan != 1 and "rowspan" not in already_processed:
            already_processed.add("rowspan")
            outfile.write(' rowspan="%s"' % self.gds_format_integer(self.rowspan, input_name="rowspan"))
        if self.colspan != 1 and "colspan" not in already_processed:
            already_processed.add("colspan")
            outfile.write(' colspan="%s"' % self.gds_format_integer(self.colspan, input_name="colspan"))

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeTH",
        fromsubclass_=False,
        pretty_print=True,
    ):
        super(typeTH, self).exportChildren(
            outfile, level, namespaceprefix_, namespacedef_, name_, True, pretty_print=pretty_print
        )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        if node.text is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", node.text)
            self.content_.append(obj_)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("rowspan", node)
        if value is not None and "rowspan" not in already_processed:
            already_processed.add("rowspan")
            self.rowspan = self.gds_parse_integer(value, node, "rowspan")
            if self.rowspan < 0:
                raise_parse_error(node, "Invalid NonNegativeInteger")
            self.validate_typeNumber(self.rowspan)  # validate type typeNumber
        value = find_attr_value_("colspan", node)
        if value is not None and "colspan" not in already_processed:
            already_processed.add("colspan")
            self.colspan = self.gds_parse_integer(value, node, "colspan")
            if self.colspan < 0:
                raise_parse_error(node, "Invalid NonNegativeInteger")
            self.validate_typeNumber(self.colspan)  # validate type typeNumber
        super(typeTH, self).buildAttributes(node, attrs, already_processed)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if not fromsubclass_ and child_.tail is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", child_.tail)
            self.content_.append(obj_)
        super(typeTH, self).buildChildren(child_, node, nodeName_, True)
        pass


# end class typeTH


class typeTD(typeFlow):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = typeFlow

    def __init__(
        self,
        p=None,
        ul=None,
        ol=None,
        table=None,
        a=None,
        i=None,
        b=None,
        sub=None,
        sup=None,
        rowspan="1",
        colspan="1",
        valueOf_=None,
        mixedclass_=None,
        content_=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        super(typeTD, self).__init__(p, ul, ol, table, a, i, b, sub, sup, valueOf_, mixedclass_, content_, **kwargs_)
        self.rowspan = _cast(None, rowspan)
        self.rowspan_nsprefix_ = None
        self.colspan = _cast(None, colspan)
        self.colspan_nsprefix_ = None
        self.valueOf_ = valueOf_
        if mixedclass_ is None:
            self.mixedclass_ = MixedContainer
        else:
            self.mixedclass_ = mixedclass_
        if content_ is None:
            self.content_ = []
        else:
            self.content_ = content_
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, typeTD)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if typeTD.subclass:
            return typeTD.subclass(*args_, **kwargs_)
        else:
            return typeTD(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_rowspan(self):
        return self.rowspan

    def set_rowspan(self, rowspan):
        self.rowspan = rowspan

    def get_colspan(self):
        return self.colspan

    def set_colspan(self, colspan):
        self.colspan = colspan

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def validate_typeNumber(self, value):
        # Validate type com:typeNumber, a restriction on xs:nonNegativeInteger.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, int):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (int)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            if not self.gds_validate_simple_patterns(self.validate_typeNumber_patterns_, value):
                self.gds_collector_.add_message(
                    'Value "%s" does not match xsd pattern restrictions: %s'
                    % (
                        encode_str_2_3(value),
                        self.validate_typeNumber_patterns_,
                    )
                )

    validate_typeNumber_patterns_ = [["^([0-9]+)$"]]

    def hasContent_(self):
        if (
            (1 if type(self.valueOf_) in [int, float] else self.valueOf_)
            or self.content_
            or super(typeTD, self).hasContent_()
        ):
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="com:", namespacedef_="", name_="typeTD", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("typeTD")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "typeTD":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeTD")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="typeTD", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="com:", name_="typeTD"):
        super(typeTD, self).exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="typeTD")
        if self.rowspan != 1 and "rowspan" not in already_processed:
            already_processed.add("rowspan")
            outfile.write(' rowspan="%s"' % self.gds_format_integer(self.rowspan, input_name="rowspan"))
        if self.colspan != 1 and "colspan" not in already_processed:
            already_processed.add("colspan")
            outfile.write(' colspan="%s"' % self.gds_format_integer(self.colspan, input_name="colspan"))

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="com:",
        namespacedef_="",
        name_="typeTD",
        fromsubclass_=False,
        pretty_print=True,
    ):
        super(typeTD, self).exportChildren(
            outfile, level, namespaceprefix_, namespacedef_, name_, True, pretty_print=pretty_print
        )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        if node.text is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", node.text)
            self.content_.append(obj_)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("rowspan", node)
        if value is not None and "rowspan" not in already_processed:
            already_processed.add("rowspan")
            self.rowspan = self.gds_parse_integer(value, node, "rowspan")
            if self.rowspan < 0:
                raise_parse_error(node, "Invalid NonNegativeInteger")
            self.validate_typeNumber(self.rowspan)  # validate type typeNumber
        value = find_attr_value_("colspan", node)
        if value is not None and "colspan" not in already_processed:
            already_processed.add("colspan")
            self.colspan = self.gds_parse_integer(value, node, "colspan")
            if self.colspan < 0:
                raise_parse_error(node, "Invalid NonNegativeInteger")
            self.validate_typeNumber(self.colspan)  # validate type typeNumber
        super(typeTD, self).buildAttributes(node, attrs, already_processed)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if not fromsubclass_ and child_.tail is not None:
            obj_ = self.mixedclass_(MixedContainer.CategoryText, MixedContainer.TypeNone, "", child_.tail)
            self.content_.append(obj_)
        super(typeTD, self).buildChildren(child_, node, nodeName_, True)
        pass


# end class typeTD


class DescriptionType(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        Comment=None,
        Submitter=None,
        Organization=None,
        Hold=None,
        SubmissionSoftware=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.Comment = Comment
        self.Comment_nsprefix_ = None
        self.Submitter = Submitter
        self.Submitter_nsprefix_ = None
        if Organization is None:
            self.Organization = []
        else:
            self.Organization = Organization
        self.Organization_nsprefix_ = None
        self.Hold = Hold
        self.Hold_nsprefix_ = None
        self.SubmissionSoftware = SubmissionSoftware
        self.SubmissionSoftware_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, DescriptionType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if DescriptionType.subclass:
            return DescriptionType.subclass(*args_, **kwargs_)
        else:
            return DescriptionType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Comment(self):
        return self.Comment

    def set_Comment(self, Comment):
        self.Comment = Comment

    def get_Submitter(self):
        return self.Submitter

    def set_Submitter(self, Submitter):
        self.Submitter = Submitter

    def get_Organization(self):
        return self.Organization

    def set_Organization(self, Organization):
        self.Organization = Organization

    def add_Organization(self, value):
        self.Organization.append(value)

    def insert_Organization_at(self, index, value):
        self.Organization.insert(index, value)

    def replace_Organization_at(self, index, value):
        self.Organization[index] = value

    def get_Hold(self):
        return self.Hold

    def set_Hold(self, Hold):
        self.Hold = Hold

    def get_SubmissionSoftware(self):
        return self.SubmissionSoftware

    def set_SubmissionSoftware(self, SubmissionSoftware):
        self.SubmissionSoftware = SubmissionSoftware

    def hasContent_(self):
        if (
            self.Comment is not None
            or self.Submitter is not None
            or self.Organization
            or self.Hold is not None
            or self.SubmissionSoftware is not None
        ):
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="DescriptionType", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("DescriptionType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "DescriptionType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="DescriptionType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="DescriptionType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="DescriptionType"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="DescriptionType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Comment is not None:
            namespaceprefix_ = self.Comment_nsprefix_ + ":" if (UseCapturedNS_ and self.Comment_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sComment>%s</%sComment>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Comment), input_name="Comment")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Submitter is not None:
            namespaceprefix_ = self.Submitter_nsprefix_ + ":" if (UseCapturedNS_ and self.Submitter_nsprefix_) else ""
            self.Submitter.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Submitter", pretty_print=pretty_print
            )
        for Organization_ in self.Organization:
            namespaceprefix_ = (
                self.Organization_nsprefix_ + ":" if (UseCapturedNS_ and self.Organization_nsprefix_) else ""
            )
            Organization_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Organization", pretty_print=pretty_print
            )
        if self.Hold is not None:
            namespaceprefix_ = self.Hold_nsprefix_ + ":" if (UseCapturedNS_ and self.Hold_nsprefix_) else ""
            self.Hold.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Hold", pretty_print=pretty_print
            )
        if self.SubmissionSoftware is not None:
            namespaceprefix_ = (
                self.SubmissionSoftware_nsprefix_ + ":"
                if (UseCapturedNS_ and self.SubmissionSoftware_nsprefix_)
                else ""
            )
            self.SubmissionSoftware.export(
                outfile,
                level,
                namespaceprefix_,
                namespacedef_="",
                name_="SubmissionSoftware",
                pretty_print=pretty_print,
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Comment":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Comment")
            value_ = self.gds_validate_string(value_, node, "Comment")
            self.Comment = value_
            self.Comment_nsprefix_ = child_.prefix
        elif nodeName_ == "Submitter":
            obj_ = typeAccount.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Submitter = obj_
            obj_.original_tagname_ = "Submitter"
        elif nodeName_ == "Organization":
            obj_ = typeOrganization.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Organization.append(obj_)
            obj_.original_tagname_ = "Organization"
        elif nodeName_ == "Hold":
            obj_ = HoldType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Hold = obj_
            obj_.original_tagname_ = "Hold"
        elif nodeName_ == "SubmissionSoftware":
            obj_ = SubmissionSoftwareType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.SubmissionSoftware = obj_
            obj_.original_tagname_ = "SubmissionSoftware"


# end class DescriptionType


class HoldType(GeneratedsSuper):
    """All data in this submission is requested to be publicly released on or
    after @release_date"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, release_date=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if isinstance(release_date, BaseStrType_):
            initvalue_ = datetime_.datetime.strptime(release_date, "%Y-%m-%d").date()
        else:
            initvalue_ = release_date
        self.release_date = initvalue_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, HoldType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if HoldType.subclass:
            return HoldType.subclass(*args_, **kwargs_)
        else:
            return HoldType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_release_date(self):
        return self.release_date

    def set_release_date(self, release_date):
        self.release_date = release_date

    def hasContent_(self):
        if ():
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="HoldType", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("HoldType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "HoldType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="HoldType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="HoldType", pretty_print=pretty_print
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="HoldType"):
        if self.release_date is not None and "release_date" not in already_processed:
            already_processed.add("release_date")
            outfile.write(' release_date="%s"' % self.gds_format_date(self.release_date, input_name="release_date"))

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="HoldType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("release_date", node)
        if value is not None and "release_date" not in already_processed:
            already_processed.add("release_date")
            try:
                self.release_date = self.gds_parse_date(value)
            except ValueError as exp:
                raise ValueError("Bad date attribute (release_date): %s" % exp)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class HoldType


class SubmissionSoftwareType(GeneratedsSuper):
    """Name of the third-party or in-house software, used by submitter to
    prepare this submission
    Version of the submission software"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, version=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.version = _cast(None, version)
        self.version_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, SubmissionSoftwareType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if SubmissionSoftwareType.subclass:
            return SubmissionSoftwareType.subclass(*args_, **kwargs_)
        else:
            return SubmissionSoftwareType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_version(self):
        return self.version

    def set_version(self, version):
        self.version = version

    def hasContent_(self):
        if ():
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="", namespacedef_="", name_="SubmissionSoftwareType", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("SubmissionSoftwareType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "SubmissionSoftwareType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="SubmissionSoftwareType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile,
                level + 1,
                namespaceprefix_,
                namespacedef_,
                name_="SubmissionSoftwareType",
                pretty_print=pretty_print,
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="SubmissionSoftwareType"):
        if self.version is not None and "version" not in already_processed:
            already_processed.add("version")
            outfile.write(
                " version=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.version), input_name="version")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="SubmissionSoftwareType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("version", node)
        if value is not None and "version" not in already_processed:
            already_processed.add("version")
            self.version = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class SubmissionSoftwareType


class ActionType(GeneratedsSuper):
    """Action is what to do - either process new submission (ProcessFile) or
    change status (ChangeStatus) of the existing one."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        action_id=None,
        submitter_tracking_id=None,
        AddFiles=None,
        AddData=None,
        ChangeStatus=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.action_id = _cast(None, action_id)
        self.action_id_nsprefix_ = None
        self.submitter_tracking_id = _cast(None, submitter_tracking_id)
        self.submitter_tracking_id_nsprefix_ = None
        self.AddFiles = AddFiles
        self.AddFiles_nsprefix_ = None
        self.AddData = AddData
        self.AddData_nsprefix_ = None
        self.ChangeStatus = ChangeStatus
        self.ChangeStatus_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, ActionType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if ActionType.subclass:
            return ActionType.subclass(*args_, **kwargs_)
        else:
            return ActionType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_AddFiles(self):
        return self.AddFiles

    def set_AddFiles(self, AddFiles):
        self.AddFiles = AddFiles

    def get_AddData(self):
        return self.AddData

    def set_AddData(self, AddData):
        self.AddData = AddData

    def get_ChangeStatus(self):
        return self.ChangeStatus

    def set_ChangeStatus(self, ChangeStatus):
        self.ChangeStatus = ChangeStatus

    def get_action_id(self):
        return self.action_id

    def set_action_id(self, action_id):
        self.action_id = action_id

    def get_submitter_tracking_id(self):
        return self.submitter_tracking_id

    def set_submitter_tracking_id(self, submitter_tracking_id):
        self.submitter_tracking_id = submitter_tracking_id

    def validate_submitter_tracking_idType(self, value):
        # Validate type submitter_tracking_idType, a restriction on xs:string.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, str):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (str)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            if len(value) > 255:
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s does not match xsd maxLength restriction on submitter_tracking_idType'
                    % {"value": encode_str_2_3(value), "lineno": lineno}
                )
                result = False

    def hasContent_(self):
        if self.AddFiles is not None or self.AddData is not None or self.ChangeStatus is not None:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="ActionType", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("ActionType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "ActionType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="ActionType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="ActionType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="ActionType"):
        if self.action_id is not None and "action_id" not in already_processed:
            already_processed.add("action_id")
            outfile.write(
                " action_id=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.action_id), input_name="action_id")),)
            )
        if self.submitter_tracking_id is not None and "submitter_tracking_id" not in already_processed:
            already_processed.add("submitter_tracking_id")
            outfile.write(
                " submitter_tracking_id=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(
                            quote_attrib(self.submitter_tracking_id), input_name="submitter_tracking_id"
                        )
                    ),
                )
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="ActionType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.AddFiles is not None:
            namespaceprefix_ = self.AddFiles_nsprefix_ + ":" if (UseCapturedNS_ and self.AddFiles_nsprefix_) else ""
            self.AddFiles.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="AddFiles", pretty_print=pretty_print
            )
        if self.AddData is not None:
            namespaceprefix_ = self.AddData_nsprefix_ + ":" if (UseCapturedNS_ and self.AddData_nsprefix_) else ""
            self.AddData.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="AddData", pretty_print=pretty_print
            )
        if self.ChangeStatus is not None:
            namespaceprefix_ = (
                self.ChangeStatus_nsprefix_ + ":" if (UseCapturedNS_ and self.ChangeStatus_nsprefix_) else ""
            )
            self.ChangeStatus.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="ChangeStatus", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("action_id", node)
        if value is not None and "action_id" not in already_processed:
            already_processed.add("action_id")
            self.action_id = value
            self.action_id = " ".join(self.action_id.split())
        value = find_attr_value_("submitter_tracking_id", node)
        if value is not None and "submitter_tracking_id" not in already_processed:
            already_processed.add("submitter_tracking_id")
            self.submitter_tracking_id = value
            self.validate_submitter_tracking_idType(
                self.submitter_tracking_id
            )  # validate type submitter_tracking_idType

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "AddFiles":
            obj_ = AddFilesType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.AddFiles = obj_
            obj_.original_tagname_ = "AddFiles"
        elif nodeName_ == "AddData":
            obj_ = AddDataType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.AddData = obj_
            obj_.original_tagname_ = "AddData"
        elif nodeName_ == "ChangeStatus":
            obj_ = ChangeStatusType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.ChangeStatus = obj_
            obj_.original_tagname_ = "ChangeStatus"


# end class ActionType


class AddFilesType(GeneratedsSuper):
    """Adding a group of files to the content of particular target archive in
    given context"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        target_db=None,
        target_db_context=None,
        File=None,
        Attribute=None,
        Meta=None,
        AttributeRefId=None,
        SequenceData=None,
        Publication=None,
        Status=None,
        Identifier=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.target_db = _cast(None, target_db)
        self.target_db_nsprefix_ = None
        self.target_db_context = _cast(None, target_db_context)
        self.target_db_context_nsprefix_ = None
        if File is None:
            self.File = []
        else:
            self.File = File
        self.File_nsprefix_ = None
        if Attribute is None:
            self.Attribute = []
        else:
            self.Attribute = Attribute
        self.Attribute_nsprefix_ = None
        if Meta is None:
            self.Meta = []
        else:
            self.Meta = Meta
        self.Meta_nsprefix_ = None
        if AttributeRefId is None:
            self.AttributeRefId = []
        else:
            self.AttributeRefId = AttributeRefId
        self.AttributeRefId_nsprefix_ = None
        if SequenceData is None:
            self.SequenceData = []
        else:
            self.SequenceData = SequenceData
        self.SequenceData_nsprefix_ = None
        if Publication is None:
            self.Publication = []
        else:
            self.Publication = Publication
        self.Publication_nsprefix_ = None
        self.Status = Status
        self.Status_nsprefix_ = None
        self.Identifier = Identifier
        self.Identifier_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, AddFilesType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if AddFilesType.subclass:
            return AddFilesType.subclass(*args_, **kwargs_)
        else:
            return AddFilesType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_File(self):
        return self.File

    def set_File(self, File):
        self.File = File

    def add_File(self, value):
        self.File.append(value)

    def insert_File_at(self, index, value):
        self.File.insert(index, value)

    def replace_File_at(self, index, value):
        self.File[index] = value

    def get_Attribute(self):
        return self.Attribute

    def set_Attribute(self, Attribute):
        self.Attribute = Attribute

    def add_Attribute(self, value):
        self.Attribute.append(value)

    def insert_Attribute_at(self, index, value):
        self.Attribute.insert(index, value)

    def replace_Attribute_at(self, index, value):
        self.Attribute[index] = value

    def get_Meta(self):
        return self.Meta

    def set_Meta(self, Meta):
        self.Meta = Meta

    def add_Meta(self, value):
        self.Meta.append(value)

    def insert_Meta_at(self, index, value):
        self.Meta.insert(index, value)

    def replace_Meta_at(self, index, value):
        self.Meta[index] = value

    def get_AttributeRefId(self):
        return self.AttributeRefId

    def set_AttributeRefId(self, AttributeRefId):
        self.AttributeRefId = AttributeRefId

    def add_AttributeRefId(self, value):
        self.AttributeRefId.append(value)

    def insert_AttributeRefId_at(self, index, value):
        self.AttributeRefId.insert(index, value)

    def replace_AttributeRefId_at(self, index, value):
        self.AttributeRefId[index] = value

    def get_SequenceData(self):
        return self.SequenceData

    def set_SequenceData(self, SequenceData):
        self.SequenceData = SequenceData

    def add_SequenceData(self, value):
        self.SequenceData.append(value)

    def insert_SequenceData_at(self, index, value):
        self.SequenceData.insert(index, value)

    def replace_SequenceData_at(self, index, value):
        self.SequenceData[index] = value

    def get_Publication(self):
        return self.Publication

    def set_Publication(self, Publication):
        self.Publication = Publication

    def add_Publication(self, value):
        self.Publication.append(value)

    def insert_Publication_at(self, index, value):
        self.Publication.insert(index, value)

    def replace_Publication_at(self, index, value):
        self.Publication[index] = value

    def get_Status(self):
        return self.Status

    def set_Status(self, Status):
        self.Status = Status

    def get_Identifier(self):
        return self.Identifier

    def set_Identifier(self, Identifier):
        self.Identifier = Identifier

    def get_target_db(self):
        return self.target_db

    def set_target_db(self, target_db):
        self.target_db = target_db

    def get_target_db_context(self):
        return self.target_db_context

    def set_target_db_context(self, target_db_context):
        self.target_db_context = target_db_context

    def validate_typeTargetDb(self, value):
        # Validate type typeTargetDb, a restriction on xs:string.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, str):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (str)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            value = value
            enumerations = [
                "BioProject",
                "BioSample",
                "clinvar",
                "dbGaP",
                "WGS",
                "variation",
                "variation_submission",
                "GTR",
                "TSA",
                "CompleteGenomes",
                "dbVar",
                "SRA",
                "SRA.experiment",
                "SRA.run",
                "SP",
                "PGAP",
                "GenBank",
                "SupFiles",
                "EST",
                "GSS",
                "TPA",
            ]
            if value not in enumerations:
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s does not match xsd enumeration restriction on typeTargetDb'
                    % {"value": encode_str_2_3(value), "lineno": lineno}
                )
                result = False

    def hasContent_(self):
        if (
            self.File
            or self.Attribute
            or self.Meta
            or self.AttributeRefId
            or self.SequenceData
            or self.Publication
            or self.Status is not None
            or self.Identifier is not None
        ):
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="AddFilesType",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("AddFilesType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "AddFilesType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="AddFilesType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="AddFilesType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="AddFilesType"):
        if self.target_db is not None and "target_db" not in already_processed:
            already_processed.add("target_db")
            outfile.write(
                " target_db=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.target_db), input_name="target_db")),)
            )
        if self.target_db_context is not None and "target_db_context" not in already_processed:
            already_processed.add("target_db_context")
            outfile.write(
                " target_db_context=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(quote_attrib(self.target_db_context), input_name="target_db_context")
                    ),
                )
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="AddFilesType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        for File_ in self.File:
            namespaceprefix_ = self.File_nsprefix_ + ":" if (UseCapturedNS_ and self.File_nsprefix_) else ""
            File_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="File", pretty_print=pretty_print)
        for Attribute_ in self.Attribute:
            namespaceprefix_ = self.Attribute_nsprefix_ + ":" if (UseCapturedNS_ and self.Attribute_nsprefix_) else ""
            Attribute_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Attribute", pretty_print=pretty_print
            )
        for Meta_ in self.Meta:
            namespaceprefix_ = self.Meta_nsprefix_ + ":" if (UseCapturedNS_ and self.Meta_nsprefix_) else ""
            Meta_.export(outfile, level, namespaceprefix_, namespacedef_="", name_="Meta", pretty_print=pretty_print)
        for AttributeRefId_ in self.AttributeRefId:
            namespaceprefix_ = (
                self.AttributeRefId_nsprefix_ + ":" if (UseCapturedNS_ and self.AttributeRefId_nsprefix_) else ""
            )
            AttributeRefId_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="AttributeRefId", pretty_print=pretty_print
            )
        for SequenceData_ in self.SequenceData:
            namespaceprefix_ = (
                self.SequenceData_nsprefix_ + ":" if (UseCapturedNS_ and self.SequenceData_nsprefix_) else ""
            )
            SequenceData_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="SequenceData", pretty_print=pretty_print
            )
        for Publication_ in self.Publication:
            namespaceprefix_ = (
                self.Publication_nsprefix_ + ":" if (UseCapturedNS_ and self.Publication_nsprefix_) else ""
            )
            Publication_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Publication", pretty_print=pretty_print
            )
        if self.Status is not None:
            namespaceprefix_ = self.Status_nsprefix_ + ":" if (UseCapturedNS_ and self.Status_nsprefix_) else ""
            self.Status.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Status", pretty_print=pretty_print
            )
        if self.Identifier is not None:
            namespaceprefix_ = self.Identifier_nsprefix_ + ":" if (UseCapturedNS_ and self.Identifier_nsprefix_) else ""
            self.Identifier.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Identifier", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("target_db", node)
        if value is not None and "target_db" not in already_processed:
            already_processed.add("target_db")
            self.target_db = value
            self.validate_typeTargetDb(self.target_db)  # validate type typeTargetDb
        value = find_attr_value_("target_db_context", node)
        if value is not None and "target_db_context" not in already_processed:
            already_processed.add("target_db_context")
            self.target_db_context = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "File":
            obj_ = FileType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.File.append(obj_)
            obj_.original_tagname_ = "File"
        elif nodeName_ == "Attribute":
            obj_ = typeFileAttribute.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Attribute.append(obj_)
            obj_.original_tagname_ = "Attribute"
        elif nodeName_ == "Meta":
            class_obj_ = self.get_class_obj_(child_, typeInlineData)
            obj_ = class_obj_.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Meta.append(obj_)
            obj_.original_tagname_ = "Meta"
        elif nodeName_ == "AttributeRefId":
            obj_ = typeFileAttributeRefId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.AttributeRefId.append(obj_)
            obj_.original_tagname_ = "AttributeRefId"
        elif nodeName_ == "SequenceData":
            obj_ = typeSequenceData.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.SequenceData.append(obj_)
            obj_.original_tagname_ = "SequenceData"
        elif nodeName_ == "Publication":
            obj_ = typePublication.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Publication.append(obj_)
            obj_.original_tagname_ = "Publication"
        elif nodeName_ == "Status":
            obj_ = typeReleaseStatus.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Status = obj_
            obj_.original_tagname_ = "Status"
        elif nodeName_ == "Identifier":
            obj_ = typeIdentifier.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Identifier = obj_
            obj_.original_tagname_ = "Identifier"


# end class AddFilesType


class FileType(typeFile):
    """File label the use of which is specific to target database. For example,
    for dbGaP genotype files, it can represent anonymized ids of dbGaP
    samples"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = typeFile

    def __init__(
        self,
        file_path=None,
        file_id=None,
        cloud_id=None,
        md5=None,
        crc32=None,
        content_type=None,
        DataType=None,
        target_db_label=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        super(FileType, self).__init__(file_path, file_id, cloud_id, md5, crc32, content_type, DataType, **kwargs_)
        self.target_db_label = _cast(None, target_db_label)
        self.target_db_label_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, FileType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if FileType.subclass:
            return FileType.subclass(*args_, **kwargs_)
        else:
            return FileType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_target_db_label(self):
        return self.target_db_label

    def set_target_db_label(self, target_db_label):
        self.target_db_label = target_db_label

    def hasContent_(self):
        if super(FileType, self).hasContent_():
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="FileType", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("FileType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "FileType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="FileType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="FileType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="FileType"):
        super(FileType, self).exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="FileType")
        if self.target_db_label is not None and "target_db_label" not in already_processed:
            already_processed.add("target_db_label")
            outfile.write(
                " target_db_label=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(quote_attrib(self.target_db_label), input_name="target_db_label")
                    ),
                )
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="FileType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        super(FileType, self).exportChildren(
            outfile, level, namespaceprefix_, namespacedef_, name_, True, pretty_print=pretty_print
        )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("target_db_label", node)
        if value is not None and "target_db_label" not in already_processed:
            already_processed.add("target_db_label")
            self.target_db_label = value
        super(FileType, self).buildAttributes(node, attrs, already_processed)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        super(FileType, self).buildChildren(child_, node, nodeName_, True)
        pass


# end class FileType


class AddDataType(GeneratedsSuper):
    """Adding a group of data objects, inlined into the submission"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        target_db=None,
        target_db_context=None,
        Data=None,
        Attribute=None,
        AttributeRefId=None,
        SequenceData=None,
        Publication=None,
        Status=None,
        Identifier=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.target_db = _cast(None, target_db)
        self.target_db_nsprefix_ = None
        self.target_db_context = _cast(None, target_db_context)
        self.target_db_context_nsprefix_ = None
        self.Data = Data
        self.Data_nsprefix_ = None
        if Attribute is None:
            self.Attribute = []
        else:
            self.Attribute = Attribute
        self.Attribute_nsprefix_ = None
        if AttributeRefId is None:
            self.AttributeRefId = []
        else:
            self.AttributeRefId = AttributeRefId
        self.AttributeRefId_nsprefix_ = None
        if SequenceData is None:
            self.SequenceData = []
        else:
            self.SequenceData = SequenceData
        self.SequenceData_nsprefix_ = None
        if Publication is None:
            self.Publication = []
        else:
            self.Publication = Publication
        self.Publication_nsprefix_ = None
        self.Status = Status
        self.Status_nsprefix_ = None
        self.Identifier = Identifier
        self.Identifier_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, AddDataType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if AddDataType.subclass:
            return AddDataType.subclass(*args_, **kwargs_)
        else:
            return AddDataType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Data(self):
        return self.Data

    def set_Data(self, Data):
        self.Data = Data

    def get_Attribute(self):
        return self.Attribute

    def set_Attribute(self, Attribute):
        self.Attribute = Attribute

    def add_Attribute(self, value):
        self.Attribute.append(value)

    def insert_Attribute_at(self, index, value):
        self.Attribute.insert(index, value)

    def replace_Attribute_at(self, index, value):
        self.Attribute[index] = value

    def get_AttributeRefId(self):
        return self.AttributeRefId

    def set_AttributeRefId(self, AttributeRefId):
        self.AttributeRefId = AttributeRefId

    def add_AttributeRefId(self, value):
        self.AttributeRefId.append(value)

    def insert_AttributeRefId_at(self, index, value):
        self.AttributeRefId.insert(index, value)

    def replace_AttributeRefId_at(self, index, value):
        self.AttributeRefId[index] = value

    def get_SequenceData(self):
        return self.SequenceData

    def set_SequenceData(self, SequenceData):
        self.SequenceData = SequenceData

    def add_SequenceData(self, value):
        self.SequenceData.append(value)

    def insert_SequenceData_at(self, index, value):
        self.SequenceData.insert(index, value)

    def replace_SequenceData_at(self, index, value):
        self.SequenceData[index] = value

    def get_Publication(self):
        return self.Publication

    def set_Publication(self, Publication):
        self.Publication = Publication

    def add_Publication(self, value):
        self.Publication.append(value)

    def insert_Publication_at(self, index, value):
        self.Publication.insert(index, value)

    def replace_Publication_at(self, index, value):
        self.Publication[index] = value

    def get_Status(self):
        return self.Status

    def set_Status(self, Status):
        self.Status = Status

    def get_Identifier(self):
        return self.Identifier

    def set_Identifier(self, Identifier):
        self.Identifier = Identifier

    def get_target_db(self):
        return self.target_db

    def set_target_db(self, target_db):
        self.target_db = target_db

    def get_target_db_context(self):
        return self.target_db_context

    def set_target_db_context(self, target_db_context):
        self.target_db_context = target_db_context

    def validate_typeTargetDb(self, value):
        # Validate type typeTargetDb, a restriction on xs:string.
        if value is not None and Validate_simpletypes_ and self.gds_collector_ is not None:
            if not isinstance(value, str):
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s is not of the correct base simple type (str)'
                    % {
                        "value": value,
                        "lineno": lineno,
                    }
                )
                return False
            value = value
            enumerations = [
                "BioProject",
                "BioSample",
                "clinvar",
                "dbGaP",
                "WGS",
                "variation",
                "variation_submission",
                "GTR",
                "TSA",
                "CompleteGenomes",
                "dbVar",
                "SRA",
                "SRA.experiment",
                "SRA.run",
                "SP",
                "PGAP",
                "GenBank",
                "SupFiles",
                "EST",
                "GSS",
                "TPA",
            ]
            if value not in enumerations:
                lineno = self.gds_get_node_lineno_()
                self.gds_collector_.add_message(
                    'Value "%(value)s"%(lineno)s does not match xsd enumeration restriction on typeTargetDb'
                    % {"value": encode_str_2_3(value), "lineno": lineno}
                )
                result = False

    def hasContent_(self):
        if (
            self.Data is not None
            or self.Attribute
            or self.AttributeRefId
            or self.SequenceData
            or self.Publication
            or self.Status is not None
            or self.Identifier is not None
        ):
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="AddDataType",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("AddDataType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "AddDataType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="AddDataType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="AddDataType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="AddDataType"):
        if self.target_db is not None and "target_db" not in already_processed:
            already_processed.add("target_db")
            outfile.write(
                " target_db=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.target_db), input_name="target_db")),)
            )
        if self.target_db_context is not None and "target_db_context" not in already_processed:
            already_processed.add("target_db_context")
            outfile.write(
                " target_db_context=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(quote_attrib(self.target_db_context), input_name="target_db_context")
                    ),
                )
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="AddDataType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Data is not None:
            namespaceprefix_ = self.Data_nsprefix_ + ":" if (UseCapturedNS_ and self.Data_nsprefix_) else ""
            self.Data.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Data", pretty_print=pretty_print
            )
        for Attribute_ in self.Attribute:
            namespaceprefix_ = self.Attribute_nsprefix_ + ":" if (UseCapturedNS_ and self.Attribute_nsprefix_) else ""
            Attribute_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Attribute", pretty_print=pretty_print
            )
        for AttributeRefId_ in self.AttributeRefId:
            namespaceprefix_ = (
                self.AttributeRefId_nsprefix_ + ":" if (UseCapturedNS_ and self.AttributeRefId_nsprefix_) else ""
            )
            AttributeRefId_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="AttributeRefId", pretty_print=pretty_print
            )
        for SequenceData_ in self.SequenceData:
            namespaceprefix_ = (
                self.SequenceData_nsprefix_ + ":" if (UseCapturedNS_ and self.SequenceData_nsprefix_) else ""
            )
            SequenceData_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="SequenceData", pretty_print=pretty_print
            )
        for Publication_ in self.Publication:
            namespaceprefix_ = (
                self.Publication_nsprefix_ + ":" if (UseCapturedNS_ and self.Publication_nsprefix_) else ""
            )
            Publication_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Publication", pretty_print=pretty_print
            )
        if self.Status is not None:
            namespaceprefix_ = self.Status_nsprefix_ + ":" if (UseCapturedNS_ and self.Status_nsprefix_) else ""
            self.Status.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Status", pretty_print=pretty_print
            )
        if self.Identifier is not None:
            namespaceprefix_ = self.Identifier_nsprefix_ + ":" if (UseCapturedNS_ and self.Identifier_nsprefix_) else ""
            self.Identifier.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Identifier", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("target_db", node)
        if value is not None and "target_db" not in already_processed:
            already_processed.add("target_db")
            self.target_db = value
            self.validate_typeTargetDb(self.target_db)  # validate type typeTargetDb
        value = find_attr_value_("target_db_context", node)
        if value is not None and "target_db_context" not in already_processed:
            already_processed.add("target_db_context")
            self.target_db_context = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Data":
            obj_ = DataType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Data = obj_
            obj_.original_tagname_ = "Data"
        elif nodeName_ == "Attribute":
            obj_ = typeFileAttribute.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Attribute.append(obj_)
            obj_.original_tagname_ = "Attribute"
        elif nodeName_ == "AttributeRefId":
            obj_ = typeFileAttributeRefId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.AttributeRefId.append(obj_)
            obj_.original_tagname_ = "AttributeRefId"
        elif nodeName_ == "SequenceData":
            obj_ = typeSequenceData.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.SequenceData.append(obj_)
            obj_.original_tagname_ = "SequenceData"
        elif nodeName_ == "Publication":
            obj_ = typePublication.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Publication.append(obj_)
            obj_.original_tagname_ = "Publication"
        elif nodeName_ == "Status":
            obj_ = typeReleaseStatus.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Status = obj_
            obj_.original_tagname_ = "Status"
        elif nodeName_ == "Identifier":
            obj_ = typeIdentifier.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Identifier = obj_
            obj_.original_tagname_ = "Identifier"


# end class AddDataType


class DataType(typeInlineData):
    """Data label the use of which is specific to target database. Same as for
    AddFiles, but this is inlined data."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = typeInlineData

    def __init__(
        self,
        name=None,
        data_model=None,
        content_type=None,
        content_encoding=None,
        XmlContent=None,
        DataContent=None,
        target_db_label=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        super(DataType, self).__init__(
            name, data_model, content_type, content_encoding, XmlContent, DataContent, **kwargs_
        )
        self.target_db_label = _cast(None, target_db_label)
        self.target_db_label_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, DataType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if DataType.subclass:
            return DataType.subclass(*args_, **kwargs_)
        else:
            return DataType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_target_db_label(self):
        return self.target_db_label

    def set_target_db_label(self, target_db_label):
        self.target_db_label = target_db_label

    def hasContent_(self):
        if super(DataType, self).hasContent_():
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="DataType", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("DataType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "DataType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="DataType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="DataType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="DataType"):
        super(DataType, self).exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="DataType")
        if self.target_db_label is not None and "target_db_label" not in already_processed:
            already_processed.add("target_db_label")
            outfile.write(
                " target_db_label=%s"
                % (
                    self.gds_encode(
                        self.gds_format_string(quote_attrib(self.target_db_label), input_name="target_db_label")
                    ),
                )
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="DataType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        super(DataType, self).exportChildren(
            outfile, level, namespaceprefix_, namespacedef_, name_, True, pretty_print=pretty_print
        )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("target_db_label", node)
        if value is not None and "target_db_label" not in already_processed:
            already_processed.add("target_db_label")
            self.target_db_label = value
        super(DataType, self).buildAttributes(node, attrs, already_processed)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        super(DataType, self).buildChildren(child_, node, nodeName_, True)
        pass


# end class DataType


class ChangeStatusType(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        Target=None,
        Release=None,
        SetReleaseDate=None,
        Suppress=None,
        Withdraw=None,
        AddComment=None,
        Identifier=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.Target = Target
        self.Target_nsprefix_ = None
        self.Release = Release
        self.Release_nsprefix_ = None
        self.SetReleaseDate = SetReleaseDate
        self.SetReleaseDate_nsprefix_ = None
        self.Suppress = Suppress
        self.Suppress_nsprefix_ = None
        self.Withdraw = Withdraw
        self.Withdraw_nsprefix_ = None
        self.AddComment = AddComment
        self.AddComment_nsprefix_ = None
        self.Identifier = Identifier
        self.Identifier_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, ChangeStatusType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if ChangeStatusType.subclass:
            return ChangeStatusType.subclass(*args_, **kwargs_)
        else:
            return ChangeStatusType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Target(self):
        return self.Target

    def set_Target(self, Target):
        self.Target = Target

    def get_Release(self):
        return self.Release

    def set_Release(self, Release):
        self.Release = Release

    def get_SetReleaseDate(self):
        return self.SetReleaseDate

    def set_SetReleaseDate(self, SetReleaseDate):
        self.SetReleaseDate = SetReleaseDate

    def get_Suppress(self):
        return self.Suppress

    def set_Suppress(self, Suppress):
        self.Suppress = Suppress

    def get_Withdraw(self):
        return self.Withdraw

    def set_Withdraw(self, Withdraw):
        self.Withdraw = Withdraw

    def get_AddComment(self):
        return self.AddComment

    def set_AddComment(self, AddComment):
        self.AddComment = AddComment

    def get_Identifier(self):
        return self.Identifier

    def set_Identifier(self, Identifier):
        self.Identifier = Identifier

    def hasContent_(self):
        if (
            self.Target is not None
            or self.Release is not None
            or self.SetReleaseDate is not None
            or self.Suppress is not None
            or self.Withdraw is not None
            or self.AddComment is not None
            or self.Identifier is not None
        ):
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="ChangeStatusType",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("ChangeStatusType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "ChangeStatusType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="ChangeStatusType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="ChangeStatusType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="ChangeStatusType"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="ChangeStatusType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Target is not None:
            namespaceprefix_ = self.Target_nsprefix_ + ":" if (UseCapturedNS_ and self.Target_nsprefix_) else ""
            self.Target.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Target", pretty_print=pretty_print
            )
        if self.Release is not None:
            namespaceprefix_ = self.Release_nsprefix_ + ":" if (UseCapturedNS_ and self.Release_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sRelease>%s</%sRelease>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Release), input_name="Release")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.SetReleaseDate is not None:
            namespaceprefix_ = (
                self.SetReleaseDate_nsprefix_ + ":" if (UseCapturedNS_ and self.SetReleaseDate_nsprefix_) else ""
            )
            self.SetReleaseDate.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="SetReleaseDate", pretty_print=pretty_print
            )
        if self.Suppress is not None:
            namespaceprefix_ = self.Suppress_nsprefix_ + ":" if (UseCapturedNS_ and self.Suppress_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sSuppress>%s</%sSuppress>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Suppress), input_name="Suppress")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Withdraw is not None:
            namespaceprefix_ = self.Withdraw_nsprefix_ + ":" if (UseCapturedNS_ and self.Withdraw_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sWithdraw>%s</%sWithdraw>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Withdraw), input_name="Withdraw")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.AddComment is not None:
            namespaceprefix_ = self.AddComment_nsprefix_ + ":" if (UseCapturedNS_ and self.AddComment_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sAddComment>%s</%sAddComment>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.AddComment), input_name="AddComment")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Identifier is not None:
            namespaceprefix_ = self.Identifier_nsprefix_ + ":" if (UseCapturedNS_ and self.Identifier_nsprefix_) else ""
            self.Identifier.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Identifier", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Target":
            obj_ = typeRefId.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Target = obj_
            obj_.original_tagname_ = "Target"
        elif nodeName_ == "Release":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Release")
            value_ = self.gds_validate_string(value_, node, "Release")
            self.Release = value_
            self.Release_nsprefix_ = child_.prefix
        elif nodeName_ == "SetReleaseDate":
            obj_ = SetReleaseDateType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.SetReleaseDate = obj_
            obj_.original_tagname_ = "SetReleaseDate"
        elif nodeName_ == "Suppress":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Suppress")
            value_ = self.gds_validate_string(value_, node, "Suppress")
            self.Suppress = value_
            self.Suppress_nsprefix_ = child_.prefix
        elif nodeName_ == "Withdraw":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Withdraw")
            value_ = self.gds_validate_string(value_, node, "Withdraw")
            self.Withdraw = value_
            self.Withdraw_nsprefix_ = child_.prefix
        elif nodeName_ == "AddComment":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "AddComment")
            value_ = self.gds_validate_string(value_, node, "AddComment")
            self.AddComment = value_
            self.AddComment_nsprefix_ = child_.prefix
        elif nodeName_ == "Identifier":
            obj_ = typeIdentifier.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Identifier = obj_
            obj_.original_tagname_ = "Identifier"


# end class ChangeStatusType


class SetReleaseDateType(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, release_date=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if isinstance(release_date, BaseStrType_):
            initvalue_ = datetime_.datetime.strptime(release_date, "%Y-%m-%d").date()
        else:
            initvalue_ = release_date
        self.release_date = initvalue_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, SetReleaseDateType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if SetReleaseDateType.subclass:
            return SetReleaseDateType.subclass(*args_, **kwargs_)
        else:
            return SetReleaseDateType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_release_date(self):
        return self.release_date

    def set_release_date(self, release_date):
        self.release_date = release_date

    def hasContent_(self):
        if ():
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="", namespacedef_="", name_="SetReleaseDateType", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("SetReleaseDateType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "SetReleaseDateType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="SetReleaseDateType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile,
                level + 1,
                namespaceprefix_,
                namespacedef_,
                name_="SetReleaseDateType",
                pretty_print=pretty_print,
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="SetReleaseDateType"):
        if self.release_date is not None and "release_date" not in already_processed:
            already_processed.add("release_date")
            outfile.write(' release_date="%s"' % self.gds_format_date(self.release_date, input_name="release_date"))

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="SetReleaseDateType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("release_date", node)
        if value is not None and "release_date" not in already_processed:
            already_processed.add("release_date")
            try:
                self.release_date = self.gds_parse_date(value)
            except ValueError as exp:
                raise ValueError("Bad date attribute (release_date): %s" % exp)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class SetReleaseDateType


class XmlContentType(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, anytypeobjs_=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if anytypeobjs_ is None:
            self.anytypeobjs_ = []
        else:
            self.anytypeobjs_ = anytypeobjs_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, XmlContentType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if XmlContentType.subclass:
            return XmlContentType.subclass(*args_, **kwargs_)
        else:
            return XmlContentType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_anytypeobjs_(self):
        return self.anytypeobjs_

    def set_anytypeobjs_(self, anytypeobjs_):
        self.anytypeobjs_ = anytypeobjs_

    def add_anytypeobjs_(self, value):
        self.anytypeobjs_.append(value)

    def insert_anytypeobjs_(self, index, value):
        self._anytypeobjs_[index] = value

    def hasContent_(self):
        if self.anytypeobjs_:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="XmlContentType", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("XmlContentType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "XmlContentType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="XmlContentType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="XmlContentType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="XmlContentType"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="XmlContentType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if not fromsubclass_:
            for obj_ in self.anytypeobjs_:
                showIndent(outfile, level, pretty_print)
                outfile.write(obj_)
                outfile.write("\n")

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        content_ = self.gds_build_any(child_, "XmlContentType")
        self.add_anytypeobjs_(content_)


# end class XmlContentType


class NameType(GeneratedsSuper):
    """Full NameName abbreviation"""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, abbr=None, valueOf_=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.abbr = _cast(None, abbr)
        self.abbr_nsprefix_ = None
        self.valueOf_ = valueOf_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, NameType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if NameType.subclass:
            return NameType.subclass(*args_, **kwargs_)
        else:
            return NameType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_abbr(self):
        return self.abbr

    def set_abbr(self, abbr):
        self.abbr = abbr

    def get_valueOf_(self):
        return self.valueOf_

    def set_valueOf_(self, valueOf_):
        self.valueOf_ = valueOf_

    def hasContent_(self):
        if 1 if type(self.valueOf_) in [int, float] else self.valueOf_:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="NameType", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("NameType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "NameType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="NameType")
        if self.hasContent_():
            outfile.write(">")
            outfile.write(self.convert_unicode(self.valueOf_))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="NameType", pretty_print=pretty_print
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="NameType"):
        if self.abbr is not None and "abbr" not in already_processed:
            already_processed.add("abbr")
            outfile.write(
                " abbr=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.abbr), input_name="abbr")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="NameType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        self.valueOf_ = get_all_text_(node)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("abbr", node)
        if value is not None and "abbr" not in already_processed:
            already_processed.add("abbr")
            self.abbr = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class NameType


class SequenceType(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, id=None, type_=None, only_one=None, Attribute=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.id = _cast(None, id)
        self.id_nsprefix_ = None
        self.type_ = _cast(None, type_)
        self.type__nsprefix_ = None
        self.only_one = _cast(None, only_one)
        self.only_one_nsprefix_ = None
        if Attribute is None:
            self.Attribute = []
        else:
            self.Attribute = Attribute
        self.Attribute_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, SequenceType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if SequenceType.subclass:
            return SequenceType.subclass(*args_, **kwargs_)
        else:
            return SequenceType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Attribute(self):
        return self.Attribute

    def set_Attribute(self, Attribute):
        self.Attribute = Attribute

    def add_Attribute(self, value):
        self.Attribute.append(value)

    def insert_Attribute_at(self, index, value):
        self.Attribute.insert(index, value)

    def replace_Attribute_at(self, index, value):
        self.Attribute[index] = value

    def get_id(self):
        return self.id

    def set_id(self, id):
        self.id = id

    def get_type(self):
        return self.type_

    def set_type(self, type_):
        self.type_ = type_

    def get_only_one(self):
        return self.only_one

    def set_only_one(self, only_one):
        self.only_one = only_one

    def hasContent_(self):
        if self.Attribute:
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="SequenceType", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("SequenceType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "SequenceType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="SequenceType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="SequenceType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="SequenceType"):
        if self.id is not None and "id" not in already_processed:
            already_processed.add("id")
            outfile.write(" id=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.id), input_name="id")),))
        if self.type_ is not None and "type_" not in already_processed:
            already_processed.add("type_")
            outfile.write(
                " type=%s" % (self.gds_encode(self.gds_format_string(quote_attrib(self.type_), input_name="type")),)
            )
        if self.only_one is not None and "only_one" not in already_processed:
            already_processed.add("only_one")
            outfile.write(
                " only_one=%s"
                % (self.gds_encode(self.gds_format_string(quote_attrib(self.only_one), input_name="only_one")),)
            )

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="SequenceType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        for Attribute_ in self.Attribute:
            namespaceprefix_ = self.Attribute_nsprefix_ + ":" if (UseCapturedNS_ and self.Attribute_nsprefix_) else ""
            Attribute_.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Attribute", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("id", node)
        if value is not None and "id" not in already_processed:
            already_processed.add("id")
            self.id = value
        value = find_attr_value_("type", node)
        if value is not None and "type" not in already_processed:
            already_processed.add("type")
            self.type_ = value
        value = find_attr_value_("only_one", node)
        if value is not None and "only_one" not in already_processed:
            already_processed.add("only_one")
            self.only_one = value

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Attribute":
            obj_ = typeFileAttribute.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Attribute.append(obj_)
            obj_.original_tagname_ = "Attribute"


# end class SequenceType


class SetReleaseDateType1(GeneratedsSuper):
    """Release on or after specific date."""

    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, release_date=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        if isinstance(release_date, BaseStrType_):
            initvalue_ = datetime_.datetime.strptime(release_date, "%Y-%m-%d").date()
        else:
            initvalue_ = release_date
        self.release_date = initvalue_

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, SetReleaseDateType1)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if SetReleaseDateType1.subclass:
            return SetReleaseDateType1.subclass(*args_, **kwargs_)
        else:
            return SetReleaseDateType1(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_release_date(self):
        return self.release_date

    def set_release_date(self, release_date):
        self.release_date = release_date

    def hasContent_(self):
        if ():
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="", namespacedef_="", name_="SetReleaseDateType1", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("SetReleaseDateType1")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "SetReleaseDateType1":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="SetReleaseDateType1")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile,
                level + 1,
                namespaceprefix_,
                namespacedef_,
                name_="SetReleaseDateType1",
                pretty_print=pretty_print,
            )
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="SetReleaseDateType1"):
        if self.release_date is not None and "release_date" not in already_processed:
            already_processed.add("release_date")
            outfile.write(' release_date="%s"' % self.gds_format_date(self.release_date, input_name="release_date"))

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="SetReleaseDateType1",
        fromsubclass_=False,
        pretty_print=True,
    ):
        pass

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        value = find_attr_value_("release_date", node)
        if value is not None and "release_date" not in already_processed:
            already_processed.add("release_date")
            try:
                self.release_date = self.gds_parse_date(value)
            except ValueError as exp:
                raise ValueError("Bad date attribute (release_date): %s" % exp)

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        pass


# end class SetReleaseDateType1


class AuthorType(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, Name=None, Consortium=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.Name = Name
        self.Name_nsprefix_ = None
        self.Consortium = Consortium
        self.Consortium_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, AuthorType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if AuthorType.subclass:
            return AuthorType.subclass(*args_, **kwargs_)
        else:
            return AuthorType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Name(self):
        return self.Name

    def set_Name(self, Name):
        self.Name = Name

    def get_Consortium(self):
        return self.Consortium

    def set_Consortium(self, Consortium):
        self.Consortium = Consortium

    def hasContent_(self):
        if self.Name is not None or self.Consortium is not None:
            return True
        else:
            return False

    def export(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="AuthorType",
        pretty_print=True,
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("AuthorType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "AuthorType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="AuthorType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="AuthorType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="AuthorType"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_=' xmlns:com="SP.common" ',
        name_="AuthorType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Name is not None:
            namespaceprefix_ = self.Name_nsprefix_ + ":" if (UseCapturedNS_ and self.Name_nsprefix_) else ""
            self.Name.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Name", pretty_print=pretty_print
            )
        if self.Consortium is not None:
            namespaceprefix_ = self.Consortium_nsprefix_ + ":" if (UseCapturedNS_ and self.Consortium_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sConsortium>%s</%sConsortium>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Consortium), input_name="Consortium")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Name":
            obj_ = typeAuthorName.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Name = obj_
            obj_.original_tagname_ = "Name"
        elif nodeName_ == "Consortium":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Consortium")
            value_ = self.gds_validate_string(value_, node, "Consortium")
            self.Consortium = value_
            self.Consortium_nsprefix_ = child_.prefix


# end class AuthorType


class StructuredCitationType(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(self, Title=None, Journal=None, gds_collector_=None, **kwargs_):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.Title = Title
        self.Title_nsprefix_ = None
        self.Journal = Journal
        self.Journal_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, StructuredCitationType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if StructuredCitationType.subclass:
            return StructuredCitationType.subclass(*args_, **kwargs_)
        else:
            return StructuredCitationType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_Title(self):
        return self.Title

    def set_Title(self, Title):
        self.Title = Title

    def get_Journal(self):
        return self.Journal

    def set_Journal(self, Journal):
        self.Journal = Journal

    def hasContent_(self):
        if self.Title is not None or self.Journal is not None:
            return True
        else:
            return False

    def export(
        self, outfile, level, namespaceprefix_="", namespacedef_="", name_="StructuredCitationType", pretty_print=True
    ):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("StructuredCitationType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "StructuredCitationType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="StructuredCitationType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile,
                level + 1,
                namespaceprefix_,
                namespacedef_,
                name_="StructuredCitationType",
                pretty_print=pretty_print,
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="StructuredCitationType"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="StructuredCitationType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.Title is not None:
            namespaceprefix_ = self.Title_nsprefix_ + ":" if (UseCapturedNS_ and self.Title_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sTitle>%s</%sTitle>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Title), input_name="Title")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Journal is not None:
            namespaceprefix_ = self.Journal_nsprefix_ + ":" if (UseCapturedNS_ and self.Journal_nsprefix_) else ""
            self.Journal.export(
                outfile, level, namespaceprefix_, namespacedef_="", name_="Journal", pretty_print=pretty_print
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "Title":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Title")
            value_ = self.gds_validate_string(value_, node, "Title")
            self.Title = value_
            self.Title_nsprefix_ = child_.prefix
        elif nodeName_ == "Journal":
            obj_ = JournalType.factory(parent_object_=self)
            obj_.build(child_, gds_collector_=gds_collector_)
            self.Journal = obj_
            obj_.original_tagname_ = "Journal"


# end class StructuredCitationType


class JournalType(GeneratedsSuper):
    __hash__ = GeneratedsSuper.__hash__
    subclass = None
    superclass = None

    def __init__(
        self,
        JournalTitle=None,
        Year=None,
        Volume=None,
        Issue=None,
        PagesFrom=None,
        PagesTo=None,
        gds_collector_=None,
        **kwargs_,
    ):
        self.gds_collector_ = gds_collector_
        self.gds_elementtree_node_ = None
        self.original_tagname_ = None
        self.parent_object_ = kwargs_.get("parent_object_")
        self.ns_prefix_ = None
        self.JournalTitle = JournalTitle
        self.JournalTitle_nsprefix_ = None
        self.Year = Year
        self.Year_nsprefix_ = None
        self.Volume = Volume
        self.Volume_nsprefix_ = None
        self.Issue = Issue
        self.Issue_nsprefix_ = None
        self.PagesFrom = PagesFrom
        self.PagesFrom_nsprefix_ = None
        self.PagesTo = PagesTo
        self.PagesTo_nsprefix_ = None

    def factory(*args_, **kwargs_):
        if CurrentSubclassModule_ is not None:
            subclass = getSubclassFromModule_(CurrentSubclassModule_, JournalType)
            if subclass is not None:
                return subclass(*args_, **kwargs_)
        if JournalType.subclass:
            return JournalType.subclass(*args_, **kwargs_)
        else:
            return JournalType(*args_, **kwargs_)

    factory = staticmethod(factory)

    def get_ns_prefix_(self):
        return self.ns_prefix_

    def set_ns_prefix_(self, ns_prefix):
        self.ns_prefix_ = ns_prefix

    def get_JournalTitle(self):
        return self.JournalTitle

    def set_JournalTitle(self, JournalTitle):
        self.JournalTitle = JournalTitle

    def get_Year(self):
        return self.Year

    def set_Year(self, Year):
        self.Year = Year

    def get_Volume(self):
        return self.Volume

    def set_Volume(self, Volume):
        self.Volume = Volume

    def get_Issue(self):
        return self.Issue

    def set_Issue(self, Issue):
        self.Issue = Issue

    def get_PagesFrom(self):
        return self.PagesFrom

    def set_PagesFrom(self, PagesFrom):
        self.PagesFrom = PagesFrom

    def get_PagesTo(self):
        return self.PagesTo

    def set_PagesTo(self, PagesTo):
        self.PagesTo = PagesTo

    def hasContent_(self):
        if (
            self.JournalTitle is not None
            or self.Year is not None
            or self.Volume is not None
            or self.Issue is not None
            or self.PagesFrom is not None
            or self.PagesTo is not None
        ):
            return True
        else:
            return False

    def export(self, outfile, level, namespaceprefix_="", namespacedef_="", name_="JournalType", pretty_print=True):
        imported_ns_def_ = GenerateDSNamespaceDefs_.get("JournalType")
        if imported_ns_def_ is not None:
            namespacedef_ = imported_ns_def_
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.original_tagname_ is not None and name_ == "JournalType":
            name_ = self.original_tagname_
        if UseCapturedNS_ and self.ns_prefix_:
            namespaceprefix_ = self.ns_prefix_ + ":"
        showIndent(outfile, level, pretty_print)
        outfile.write(
            "<%s%s%s"
            % (
                namespaceprefix_,
                name_,
                namespacedef_ and " " + namespacedef_ or "",
            )
        )
        already_processed = set()
        self.exportAttributes(outfile, level, already_processed, namespaceprefix_, name_="JournalType")
        if self.hasContent_():
            outfile.write(">%s" % (eol_,))
            self.exportChildren(
                outfile, level + 1, namespaceprefix_, namespacedef_, name_="JournalType", pretty_print=pretty_print
            )
            showIndent(outfile, level, pretty_print)
            outfile.write("</%s%s>%s" % (namespaceprefix_, name_, eol_))
        else:
            outfile.write("/>%s" % (eol_,))

    def exportAttributes(self, outfile, level, already_processed, namespaceprefix_="", name_="JournalType"):
        pass

    def exportChildren(
        self,
        outfile,
        level,
        namespaceprefix_="",
        namespacedef_="",
        name_="JournalType",
        fromsubclass_=False,
        pretty_print=True,
    ):
        if pretty_print:
            eol_ = "\n"
        else:
            eol_ = ""
        if self.JournalTitle is not None:
            namespaceprefix_ = (
                self.JournalTitle_nsprefix_ + ":" if (UseCapturedNS_ and self.JournalTitle_nsprefix_) else ""
            )
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sJournalTitle>%s</%sJournalTitle>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.JournalTitle), input_name="JournalTitle")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Year is not None:
            namespaceprefix_ = self.Year_nsprefix_ + ":" if (UseCapturedNS_ and self.Year_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sYear>%s</%sYear>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Year), input_name="Year")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Volume is not None:
            namespaceprefix_ = self.Volume_nsprefix_ + ":" if (UseCapturedNS_ and self.Volume_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sVolume>%s</%sVolume>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Volume), input_name="Volume")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.Issue is not None:
            namespaceprefix_ = self.Issue_nsprefix_ + ":" if (UseCapturedNS_ and self.Issue_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sIssue>%s</%sIssue>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.Issue), input_name="Issue")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.PagesFrom is not None:
            namespaceprefix_ = self.PagesFrom_nsprefix_ + ":" if (UseCapturedNS_ and self.PagesFrom_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sPagesFrom>%s</%sPagesFrom>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.PagesFrom), input_name="PagesFrom")),
                    namespaceprefix_,
                    eol_,
                )
            )
        if self.PagesTo is not None:
            namespaceprefix_ = self.PagesTo_nsprefix_ + ":" if (UseCapturedNS_ and self.PagesTo_nsprefix_) else ""
            showIndent(outfile, level, pretty_print)
            outfile.write(
                "<%sPagesTo>%s</%sPagesTo>%s"
                % (
                    namespaceprefix_,
                    self.gds_encode(self.gds_format_string(quote_xml(self.PagesTo), input_name="PagesTo")),
                    namespaceprefix_,
                    eol_,
                )
            )

    def build(self, node, gds_collector_=None):
        self.gds_collector_ = gds_collector_
        if SaveElementTreeNode:
            self.gds_elementtree_node_ = node
        already_processed = set()
        self.ns_prefix_ = node.prefix
        self.buildAttributes(node, node.attrib, already_processed)
        for child in node:
            nodeName_ = Tag_pattern_.match(child.tag).groups()[-1]
            self.buildChildren(child, node, nodeName_, gds_collector_=gds_collector_)
        return self

    def buildAttributes(self, node, attrs, already_processed):
        pass

    def buildChildren(self, child_, node, nodeName_, fromsubclass_=False, gds_collector_=None):
        if nodeName_ == "JournalTitle":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "JournalTitle")
            value_ = self.gds_validate_string(value_, node, "JournalTitle")
            self.JournalTitle = value_
            self.JournalTitle_nsprefix_ = child_.prefix
        elif nodeName_ == "Year":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Year")
            value_ = self.gds_validate_string(value_, node, "Year")
            self.Year = value_
            self.Year_nsprefix_ = child_.prefix
        elif nodeName_ == "Volume":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Volume")
            value_ = self.gds_validate_string(value_, node, "Volume")
            self.Volume = value_
            self.Volume_nsprefix_ = child_.prefix
        elif nodeName_ == "Issue":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "Issue")
            value_ = self.gds_validate_string(value_, node, "Issue")
            self.Issue = value_
            self.Issue_nsprefix_ = child_.prefix
        elif nodeName_ == "PagesFrom":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "PagesFrom")
            value_ = self.gds_validate_string(value_, node, "PagesFrom")
            self.PagesFrom = value_
            self.PagesFrom_nsprefix_ = child_.prefix
        elif nodeName_ == "PagesTo":
            value_ = child_.text
            value_ = self.gds_parse_string(value_, node, "PagesTo")
            value_ = self.gds_validate_string(value_, node, "PagesTo")
            self.PagesTo = value_
            self.PagesTo_nsprefix_ = child_.prefix


# end class JournalType


GDSClassesMapping = {}


USAGE_TEXT = """
Usage: python <Parser>.py [ -s ] <in_xml_file>
"""


def usage():
    print(USAGE_TEXT)
    sys.exit(1)


def get_root_tag(node):
    tag = Tag_pattern_.match(node.tag).groups()[-1]
    rootClass = GDSClassesMapping.get(tag)
    if rootClass is None:
        rootClass = globals().get(tag)
    return tag, rootClass


def get_required_ns_prefix_defs(rootNode):
    """Get all name space prefix definitions required in this XML doc.
    Return a dictionary of definitions and a char string of definitions.
    """
    nsmap = {prefix: uri for node in rootNode.iter() for (prefix, uri) in node.nsmap.items() if prefix is not None}
    namespacedefs = " ".join(['xmlns:{}="{}"'.format(prefix, uri) for prefix, uri in nsmap.items()])
    return nsmap, namespacedefs


def parse(inFileName, silence=False, print_warnings=True):
    global CapturedNsmap_
    gds_collector = GdsCollector_()
    parser = None
    doc = parsexml_(inFileName, parser)
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = "Submission"
        rootClass = Submission
    rootObj = rootClass.factory()
    rootObj.build(rootNode, gds_collector_=gds_collector)
    CapturedNsmap_, namespacedefs = get_required_ns_prefix_defs(rootNode)
    if not SaveElementTreeNode:
        doc = None
        rootNode = None
    if not silence:
        sys.stdout.write('<?xml version="1.0" ?>\n')
        rootObj.export(sys.stdout, 0, name_=rootTag, namespacedef_=namespacedefs, pretty_print=True)
    if print_warnings and len(gds_collector.get_messages()) > 0:
        separator = ("-" * 50) + "\n"
        sys.stderr.write(separator)
        sys.stderr.write(
            "----- Warnings -- count: {} -----\n".format(
                len(gds_collector.get_messages()),
            )
        )
        gds_collector.write_messages(sys.stderr)
        sys.stderr.write(separator)
    return rootObj


def parseEtree(inFileName, silence=False, print_warnings=True, mapping=None, nsmap=None):
    parser = None
    doc = parsexml_(inFileName, parser)
    gds_collector = GdsCollector_()
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = "Submission"
        rootClass = Submission
    rootObj = rootClass.factory()
    rootObj.build(rootNode, gds_collector_=gds_collector)
    # Enable Python to collect the space used by the DOM.
    if mapping is None:
        mapping = {}
    rootElement = rootObj.to_etree(None, name_=rootTag, mapping_=mapping, nsmap_=nsmap)
    reverse_mapping = rootObj.gds_reverse_node_mapping(mapping)
    if not SaveElementTreeNode:
        doc = None
        rootNode = None
    if not silence:
        content = etree_.tostring(rootElement, pretty_print=True, xml_declaration=True, encoding="utf-8")
        sys.stdout.write(str(content))
        sys.stdout.write("\n")
    if print_warnings and len(gds_collector.get_messages()) > 0:
        separator = ("-" * 50) + "\n"
        sys.stderr.write(separator)
        sys.stderr.write(
            "----- Warnings -- count: {} -----\n".format(
                len(gds_collector.get_messages()),
            )
        )
        gds_collector.write_messages(sys.stderr)
        sys.stderr.write(separator)
    return rootObj, rootElement, mapping, reverse_mapping


def parseString(inString, silence=False, print_warnings=True):
    """Parse a string, create the object tree, and export it.

    Arguments:
    - inString -- A string.  This XML fragment should not start
      with an XML declaration containing an encoding.
    - silence -- A boolean.  If False, export the object.
    Returns -- The root object in the tree.
    """
    parser = None
    rootNode = parsexmlstring_(inString, parser)
    gds_collector = GdsCollector_()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = "Submission"
        rootClass = Submission
    rootObj = rootClass.factory()
    rootObj.build(rootNode, gds_collector_=gds_collector)
    if not SaveElementTreeNode:
        rootNode = None
    if not silence:
        sys.stdout.write('<?xml version="1.0" ?>\n')
        rootObj.export(sys.stdout, 0, name_=rootTag, namespacedef_="")
    if print_warnings and len(gds_collector.get_messages()) > 0:
        separator = ("-" * 50) + "\n"
        sys.stderr.write(separator)
        sys.stderr.write(
            "----- Warnings -- count: {} -----\n".format(
                len(gds_collector.get_messages()),
            )
        )
        gds_collector.write_messages(sys.stderr)
        sys.stderr.write(separator)
    return rootObj


def parseLiteral(inFileName, silence=False, print_warnings=True):
    parser = None
    doc = parsexml_(inFileName, parser)
    gds_collector = GdsCollector_()
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = "Submission"
        rootClass = Submission
    rootObj = rootClass.factory()
    rootObj.build(rootNode, gds_collector_=gds_collector)
    # Enable Python to collect the space used by the DOM.
    if not SaveElementTreeNode:
        doc = None
        rootNode = None
    if not silence:
        sys.stdout.write("#from genbank_submission import *\n\n")
        sys.stdout.write("import genbank_submission as model_\n\n")
        sys.stdout.write("rootObj = model_.rootClass(\n")
        rootObj.exportLiteral(sys.stdout, 0, name_=rootTag)
        sys.stdout.write(")\n")
    if print_warnings and len(gds_collector.get_messages()) > 0:
        separator = ("-" * 50) + "\n"
        sys.stderr.write(separator)
        sys.stderr.write(
            "----- Warnings -- count: {} -----\n".format(
                len(gds_collector.get_messages()),
            )
        )
        gds_collector.write_messages(sys.stderr)
        sys.stderr.write(separator)
    return rootObj


def main():
    args = sys.argv[1:]
    if len(args) == 1:
        parse(args[0])
    else:
        usage()


if __name__ == "__main__":
    # import pdb; pdb.set_trace()
    main()

RenameMappings_ = {}

#
# Mapping of namespaces to types defined in them
# and the file in which each is defined.
# simpleTypes are marked "ST" and complexTypes "CT".
NamespaceToDefMappings_ = {
    "SP.common": [
        ("typeNumber", "SP.common.xsd", "ST"),
        ("typeArchive", "SP.common.xsd", "ST"),
        ("typeLocalId", "SP.common.xsd", "CT"),
        ("typeSPUID", "SP.common.xsd", "CT"),
        ("typePrimaryId", "SP.common.xsd", "CT"),
        ("typeIdentifier", "SP.common.xsd", "CT"),
        ("typeRefId", "SP.common.xsd", "CT"),
        ("typeLink", "SP.common.xsd", "CT"),
        ("typeExternalLink", "SP.common.xsd", "CT"),
        ("typeAuthorSet", "SP.common.xsd", "CT"),
        ("typePublication", "SP.common.xsd", "CT"),
        ("typeDescriptor", "SP.common.xsd", "CT"),
        ("typeAddress", "SP.common.xsd", "CT"),
        ("typeName", "SP.common.xsd", "CT"),
        ("typeAuthorName", "SP.common.xsd", "CT"),
        ("typeContactInfo", "SP.common.xsd", "CT"),
        ("typeOrganism", "SP.common.xsd", "CT"),
        ("typeBlock", "SP.common.xsd", "CT"),
        ("typeInline", "SP.common.xsd", "CT"),
        ("typeFlow", "SP.common.xsd", "CT"),
        ("typeA_content", "SP.common.xsd", "CT"),
        ("typeL", "SP.common.xsd", "CT"),
        ("typeLI", "SP.common.xsd", "CT"),
        ("typeA", "SP.common.xsd", "CT"),
        ("typeTable", "SP.common.xsd", "CT"),
        ("typeCaption", "SP.common.xsd", "CT"),
        ("typeTR", "SP.common.xsd", "CT"),
        ("typeTH", "SP.common.xsd", "CT"),
        ("typeTD", "SP.common.xsd", "CT"),
    ]
}

__all__ = [
    "ActionType",
    "AddDataType",
    "AddFilesType",
    "AuthorType",
    "ChangeStatusType",
    "DataType",
    "DescriptionType",
    "FileType",
    "HoldType",
    "JournalType",
    "NameType",
    "Release",
    "SequenceType",
    "SetReleaseDateType",
    "SetReleaseDateType1",
    "StructuredCitationType",
    "Submission",
    "SubmissionSoftwareType",
    "XmlContentType",
    "typeA",
    "typeA_content",
    "typeAccount",
    "typeAddress",
    "typeAuthorName",
    "typeAuthorSet",
    "typeBlock",
    "typeCaption",
    "typeContactInfo",
    "typeDescriptor",
    "typeExternalLink",
    "typeFile",
    "typeFileAttribute",
    "typeFileAttributeRefId",
    "typeFlow",
    "typeIdentifier",
    "typeInline",
    "typeInlineData",
    "typeL",
    "typeLI",
    "typeLink",
    "typeLocalId",
    "typeName",
    "typeOrganism",
    "typeOrganization",
    "typePrimaryId",
    "typePublication",
    "typeRefId",
    "typeReleaseStatus",
    "typeSPUID",
    "typeSequenceData",
    "typeTD",
    "typeTH",
    "typeTR",
    "typeTable",
]
