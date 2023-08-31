podman build --tag github_runner . && \
podman run --security-opt=label=disable -d=true \
--env GITHUB_TOKEN=AFPDS6753IL4TZY3PPHNNZLEUWJHA \
--env REPO_URL=https://github.com/Jooorgen/madgraph4gpu \
--env GITHUB_RUNNER_TAGS=Linux,x64,a100 \
--env RUNNER_NAME=GPURunner_itscrd-a100 \
--name github_runner github_runner