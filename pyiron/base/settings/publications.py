# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

from pyiron.base.settings.generic import Settings

__author__ = "Joerg Neugebauer, Jan Janssen"
__copyright__ = (
    "Copyright 2020, Max-Planck-Institut für Eisenforschung GmbH - "
    "Computational Materials Design (CM) Department"
)
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"

s = Settings()


def list_publications(bib_format="dict"):
    """
    List the publications used in this project.

    Args:
        bib_format (str): ['dict', 'bibtex', 'apa']

    Returns:
        list: list of publications in Bibtex format.
    """

    def get_bibtex(key, value):
        total_keys = [
            "title",
            "journal",
            "volume",
            "issue",
            "number",
            "pages",
            "numpages",
            "year",
            "month",
            "publisher",
            "url",
            "doi",
            "issn",
        ]
        bibtex_str = (
                "@article{"
                + key
                + ",\n"
                + "    author={"
                + " and ".join(value["author"])
                + "},\n"
        )
        for key in total_keys:
            if key in value.keys():
                bibtex_str += "    " + key + "={" + value[key] + "},\n"
        bibtex_str += "}\n"
        return bibtex_str

    def get_apa(value):
        apa_str = " & ".join(value["author"])
        if "year" in value.keys():
            apa_str += " (" + value["year"] + "). "
        if "title" in value.keys():
            apa_str += value["title"] + ". "
        if "journal" in value.keys():
            apa_str += value["journal"] + ", "
        if "volume" in value.keys():
            apa_str += value["volume"] + ", "
        if "pages" in value.keys():
            apa_str += value["pages"] + ". "
        if "doi" in value.keys():
            apa_str += "doi: " + value["doi"] + "\n"
        return apa_str

    publication_dict = s.publication_lst
    if bib_format.lower() == "dict":
        return publication_dict
    elif bib_format.lower() == "bibtex":
        total_str = ""
        for pub in publication_dict:
            for key, value in pub.items():
                total_str += get_bibtex(key, value)
        return total_str
    elif bib_format.lower() == "apa":
        total_str = ""
        for pub in publication_dict:
            for key, value in pub.items():
                total_str += get_apa(value)
        return total_str
    else:
        raise ValueError("Supported Bibformats are ['dict', 'bibtex', 'apa']")
