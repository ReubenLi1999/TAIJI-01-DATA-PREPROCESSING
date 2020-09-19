# assume midnight at Greenwich, arbitrary date
import requests

payload = {
    'year': 2020,
    'month': 5,
    'day': 19,
    'geodetic_latitude': 51.47879,
    'geodetic_longitude': 56.93,
    'altitude': 500,  # km
    'utc': 0,  # hours
}

url = "https://atmosphere.amentum.space/jb2008"

# handle exceptions
response = requests.get(url, params=payload)
json_payload = response.json()
print(json_payload['total_mass_density']['value'])
