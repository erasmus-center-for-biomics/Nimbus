#!/bin/env python

import functools
import re
import sys
import table


@functools.total_ordering
class Variant(table.Loc):
    """
    A class to represent a VCF variant.

    This class contains the location, id, reference
    and alternate alleles in a VCF file. The id field is
    renamed to name to avoid clashes with the Python
    buildin.
    """

    __slots__ = ["id", "reference", "alternate"]

    def __init__(self, chrid, pos, name, reference, alternate):
        """Initialize a new Variant object."""
        super(Variant, self).__init__(chrid, pos - 1, pos - 1)
        self.name = name
        self.reference = reference
        self.alternate = alternate

    def __eq__(self, other):
        """Object is equal to other."""
        return (
            self.chrid,
            self.start,
            self.end,
            self.name,
            self.reference,
            self.alternate) == (
                other.chrid,
                other.start,
                other.end,
                other.name,
                other.reference,
                other.alternate)

    def __lt__(self, other):
        """Compare to other."""
        return (
            self.chrid,
            self.start,
            self.end,
            self.name,
            self.reference,
            self.alternate) < (
                other.chrid,
                other.start,
                other.end,
                other.name,
                other.reference,
                other.alternate)

    def __str__(self):
        """Represent the object as a string."""
        additional = ";id=%s;reference=%s;alternate=[%s]" % (
            self.name,
            self.reference,
            ','.join(self.alternate))
        return super(Variant, self).__str__() + additional


class Meta(object):
    """A class to represent the meta information from a VCF file."""

    # the regular expression to parse a meta line
    meta_entry_regexp = re.compile(r"^##([\w,.]+)=<?(.+)>?$")

    def __init__(self):
        self.meta = []  # a list of paired pairs to hold the information
        self.indexed = False

    # def info_entries(self):
    #     if not self.indexed:
    #         self.index()
    #     for entry in self.info:
    #         yield entry

    # def format_entries(self):
    #     if not self.indexed:
    #         self.index()
    #     for entry in self.format:
    #         yield entry

    # def index(self):
    #     self.info = []
    #     self.info_ids = []
    #     self.format = []
    #     self.format_ids = []
    #     for tag, entry in self.meta:
    #         entry = dict(entry)
    #         if tag == "INFO":
    #             self.info.append(entry)
    #             if "ID" in entry.keys():
    #                 self.info_ids.append(entry["ID"])
    #         elif tag == "FORMAT":
    #             self.format.append(entry)
    #             if "ID" in entry.keys():
    #                 self.format_ids.append(entry["ID"])
    #     self.indexed = True

    def add_from_line(self, text=""):
        """Add the content of meta line to the object."""
        self.indexed = False

        # match the header to the meta format
        text = text.rstrip()
        m = self.meta_entry_regexp.match(text)
        entry = []
        if m:
            # get the tag and content
            tag = m.group(1)
            content = m.group(2)
            content = content.rstrip(">")

            # get the pairs in the meta line
            pairs = content.split(",")
            for p in pairs:
                items = p.split("=")
                if len(items) > 2:
                    continue
                if len(items) == 2:
                    entry.append((items[0], items[1]))
                else:
                    entry.append((tag, items[0]))

            # add the meta information
            self.meta.append((tag, entry))
        else:
            raise ValueError("Could not parse meta line: %s" % text)

    def write(self, h=sys.stdout):
        for tag, content in self.meta:
            line = None

            if len(content) > 1:
                elements = []
                for a, b in content:
                    elements.append(a + '=' + b)
                content = ','.join(elements)
                line = "##%(tag)s=<%(content)s>" % {
                    "tag": tag,
                    "content": content
                }
            else:
                line = "##%(tag)s=%(content)s" % {
                    "tag": tag,
                    "content": content[0][1]
                }
            h.write("%s\n" % line)

