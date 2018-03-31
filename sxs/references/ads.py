

ads_api_url = 'https://api.adsabs.harvard.edu/v1/search/query'


def ads_bearer_token(error_if_not_found=True):
    import os.environ
    if 'ADS_TOKEN' not in os.environ:
        message = 'Environment variable "ADS_TOKEN" not found'
        if error_if_not_found:
            raise EnvironmentError(message)
        else:
            from warnings import warn
            warn(message)
            return ''
    return os.environ('ADS_TOKEN')


def authorization_header(error_if_not_found=True):
    token = ads_bearer_token(error_if_not_found)
    return 'Authorization: Bearer:{0}'.format(token)



#curl -H 'Authorization: Bearer:Co3woTdgh00kgY46RrP51sVHZBFEiml6LeQq2xNK' 'http://api.adsabs.harvard.edu/v1/search/query?q=doi:10.1111/j.1365-2966.2010.18127.x&fl=pub+volume+issue+year+page'
