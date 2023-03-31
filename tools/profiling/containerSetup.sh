# Preliminary setup
podman=${podman:-podman}
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
runnerName=GPURunner_itscrd-a100
sourceImage=nvidia/cuda:12.0.1-devel-rockylinux8
tag=githubci-cuda12.0.1-gcc11.3-clang
GitHubRunnerTags=Linux,x64,a100
githubToken=$1

# Links
runnerURL=https://github.com/actions/runner/releases/download/v2.301.1/actions-runner-linux-x64-2.301.1.tar.gz
nvidiaContainerToolkitLink=https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.repo
repoURL=https://github.com/Jooorgen/madgraph4gpu

if ! which podman > /dev/null; then
  echo "Podman not installed. Trying now ..."
  sudo yum install podman
  curl -s -L $nvidiaContainerToolkitLink > nvidia-container-runtime.repo
  sudo mv nvidia-container-runtime.repo /etc/yum-puppet.repos.d/
  sudo yum install nvidia-container-runtime

  sudo sed -i 's/^#no-cgroups = false/no-cgroups = true/;' /etc/nvidia-container-runtime/config.toml
  exit 0
fi

if $runTest; then
  # Test that container starts up
  $podman run --rm --security-opt=label=disable nvidia/cuda:11.5.0-devel-centos8 nvidia-smi || exit 1
fi

cat > entrypoint.sh << "EOF"
#!/bin/bash
RUNNER=/home/CI/actions-runner/run.sh

while true; do
  if ! pgrep -f ${RUNNER} > /dev/null 2>&1; then
    # Runner hasn't been started yet or exited because of failure / update
    ${RUNNER}
  else
    # Runner was restarted, and is running in background. Let's wait for it.
    PID=$(pgrep -f ${RUNNER}) && tail --pid=$PID -f /dev/null
  fi
  sleep 10
done
EOF

# In container:
# - install cmake, git, which
cat > containerManifest <<EOF
FROM ${sourceImage}
LABEL maintaner="Stephan/Jorgen"

# Add ARG instructions for required variables
ARG GITHUB_TOKEN
ARG RUNNER_NAME
ARG GITHUB_RUNNER_TAGS
ARG REPO_URL
ARG RUNNER_URL https://github.com/actions/runner/releases/download/v2.301.1/actions-runner-linux-x64-2.301.1.tar.gz

RUN if [ -z "\${GITHUB_TOKEN}" ]; then echo "GITHUB_TOKEN is not set" && exit 1; fi

USER root
RUN yum install -y cmake which git libicu glibc lttng-ust vim clang python39 ncurses ncurses-devel libzstd libzstd-devel findutils
RUN useradd CI
USER CI
WORKDIR /home/CI/
RUN mkdir actions-runner && cd /tmp/ && curl -o /tmp/actions-runner.tar.gz -L ${Runner_Url} && cd /home/CI/actions-runner && tar -xzf /tmp/actions-runner.tar.gz
WORKDIR /home/CI/actions-runner
RUN ./config.sh --unattended --url ${REPO_URL} --token ${GITHUB_TOKEN} --replace --name ${RUNNER_NAME} --labels ${GITHUB_RUNNER_TAGS}
COPY ./entrypoint.sh .
USER root
RUN chown CI ./entrypoint.sh && ls -l && chmod u+x ./entrypoint.sh
USER CI
CMD [ "./entrypoint.sh" ]
EOF

$podman build \
  --build-arg GITHUB_TOKEN=$githubToken \
  --build-arg RUNNER_NAME=$runnerName \
  --build-arg GITHUB_RUNNER_TAGS=$GitHubRunnerTags \
  --build-arg REPO_URL=$repoUrl \
  --build-arg RUNNER_URL=https://github.com/actions/runner/releases/download/v2.301.1/actions-runner-linux-x64-2.301.1.tar.gz \
  --tag ${tag} \
  --file containerManifest \
  || exit 1

# Run container:
# label=disable disables carrying over of SELinux labels for mounts inside the container
$podman run \
--security-opt=label=disable \
-v /cvmfs/sft.cern.ch/lcg/releases/binutils:/cvmfs/sft.cern.ch/lcg/releases/binutils:ro \
-v /cvmfs/projects.cern.ch/intelsw/oneAPI:/cvmfs/projects.cern.ch/intelsw/oneAPI:ro \
-v /cvmfs/sft.cern.ch/lcg/releases/gcc:/cvmfs/sft.cern.ch/lcg/releases/gcc:ro \
--hooks-dir=/usr/share/containers/oci/hooks.d/ \
--name github_runner ${tag}

#$podman create --security-opt=label=disable --name github_runner ${tag}
#$podman generate systemd --restart-policy=always --files -t 10 -n github_runner