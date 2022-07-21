import argparse
import requests

from oauth2client.client import GoogleCredentials

# function to get authorization bearer token for requests
def get_access_token():
    """Get access token."""

    scopes = ["https://www.googleapis.com/auth/userinfo.profile", "https://www.googleapis.com/auth/userinfo.email"]
    credentials = GoogleCredentials.get_application_default()
    credentials = credentials.create_scoped(scopes)
    return credentials.get_access_token().access_token

def call_flexible_import_entities(workspace_name, project, tsv):
    """Post entities to Terra workspace using flexibleImportEntities."""

    # rawls request URL for batchUpsert
    uri = f"https://api.firecloud.org/api/workspaces/{project}/{workspace_name}/flexibleImportEntities?async=false&deleteEmptyValues=false"
    # Get access token and and add to headers for requests.
    # -H  "accept: */*" -H  "Authorization: Bearer [token] -H "Content-Type: application/json"
    headers = {"Authorization": "Bearer " + get_access_token(), "accept": "*/*"}

	# Create file dictionary to be passed to request
    files = {'entities': open(tsv ,'rb')}

    # capture response from API and parse out status code
    response = requests.post(uri, headers=headers, files=files)
    status_code = response.status_code

    if status_code != 200:  # entities upsert fail
        print(f"ERROR: Code {status_code} returned.")
        print(response.text)
        print(response.raise_for_status())
        
    # entities upsert success
    print(f"Successfully uploaded entities." + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-w', '--workspace_name', required=True, help='name of workspace in which to make changes')
    parser.add_argument('-p', '--project', required=True, help='billing project (namespace) of workspace in which to make changes')
    parser.add_argument('-t', '--tsv', required=True, help='.tsv file formatted in load format to Terra UI')

    args = parser.parse_args()

    # call import API (firecloud)
    call_flexible_import_entities(args.workspace_name, args.project, args.tsv)