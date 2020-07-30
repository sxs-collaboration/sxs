from .api import DOI, Login

class SXSLogin(Login):
    """Use default values when creating/modifying DOI objects

    This class subclasses `Login` to add default values relevant to SXS's DOIs when creating or
    modifying DOIs on DataCite.  Specifically, we have default values for all the required values,
    except for `doi` and `url`:

      creators = [{"name": "SXS Collaboration"}]
      title = sxs.simulation_title(doi.suffix)
      publisher = "CaltechDATA"
      publicationYear = datetime.datetime.now().year
      resourceTypeGeneral = "Dataset"

    """

    prefix = "SXS PREFIX GOES HERE!!!"

    def update(self, doi):
        doi = DOI(doi)
        raise NotImplementedError("Have to add default information")

        return super().update(self, doi)
