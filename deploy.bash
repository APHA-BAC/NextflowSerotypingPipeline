docker pull jguzinski/salmonella-seq:prod
docker tag jguzinski/salmonella-seq:prod 714385292749.dkr.ecr.eu-west-1.amazonaws.com/salmonella-ec2-0-1-1:latest
aws ecr get-login-password --region eu-west-1 | docker login --username AWS --password-stdin 714385292749.dkr.ecr.eu-west-1.amazonaws.com
docker push 714385292749.dkr.ecr.eu-west-1.amazonaws.com/salmonella-ec2-0-1-1:latest
