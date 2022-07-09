import boto3 

def get_ev(file, access_key, secret_key, bucket='all-abstract-ev'):
    session = boto3.Session(
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key
    )
    s3 = session.resource("s3")
    bucket = bucket
    print("Getting", file)
    obj = s3.Object(bucket, file)
    try:
        text = obj.get()['Body'].read().decode('utf-8') 
    except Exception as e:
        print(e)
        text = "Not found"
    print("Got", file)
    return text
