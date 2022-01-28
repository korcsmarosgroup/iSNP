import os
from datetime import datetime

import docker


class DockerHelper:

    def __init__(self, navigomix_folder, display, input_folder, output_folder, image_name_with_tag, silent=True, separate=False):
        self.silent = silent
        self.image_name_with_tag = image_name_with_tag
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.docker_client_low_level = docker.APIClient(base_url=os.environ['DOCKER_HOST'])
        self.docker_client_high_level = docker.from_env()
        self.navigomix_folder = navigomix_folder
        self.display = display
        self.separate = separate
        self.container = None

    def build_docker(self, path, image_name_with_tag):
        self.display.print("  building docker image: " + image_name_with_tag)
        for output_line in self.docker_client_low_level.build(path=path, rm=True, forcerm=True, tag=image_name_with_tag, quiet=False,
                                                              nocache=False, decode=True):
            if not self.silent:
                if 'stream' in output_line:
                    self.display.print(output_line['stream'].strip())
                else:
                    self.display.print(str(output_line))  # eg. for error lines
        return image_name_with_tag

    def execute_analytical_task(self, command):
        if self.separate:
            return self.start_single_command_in_separated_container(command)
        else:
            return self.exec_command_in_long_term_container(command)

    def start_long_term_container(self):
        self.container = self.run_docker(['sleep', '864000'])  # wait for a day before timeout

    def kill_running_container(self):
        if self.container:
            self.container.remove(force=True)

    def exec_command_in_long_term_container(self, command, foreground_failure_message=''):
        if not self.silent:
            self.display.print("\n\n======= running analytical task in long-term docker container with command: %s " % command)
        (exit_code, output) = self.container.exec_run(command, environment={'PYTHONPATH': '/analytic-modules'}, stream=False, stdout=True, stderr=True, workdir='/')
        if not self.silent:
            self.display.print(str(output).replace('\\n', '\n').replace('\\\'', '\''))
        if exit_code != 0:
            if not self.silent:
                self.display.print("\n\n==== analytical task executed in long-term container exited with non-zero exit code: %d" % exit_code)
                self.display.print(foreground_failure_message)
        return exit_code

    def start_single_command_in_separated_container(self, command, run_in_foreground=True, foreground_failure_message=''):
        self.container = self.run_docker(command)

        if (run_in_foreground):
            for line in self.container.logs(stream=True):
                if not self.silent:
                    self.display.print(line.strip())
            self.container.reload()
            exit_code = self.container.attrs['State']['ExitCode']
            self.container.remove()
            self.container = None
            if exit_code != 0:
                if not self.silent:
                    self.display.print("\n\n==== docker container exited with non-zero exit code: %d" % exit_code)
                    self.display.print(foreground_failure_message)
            return exit_code

        self.container.remove()
        self.container = None
        return 0

    def run_docker(self, command):
        if not ":" in self.image_name_with_tag:
            raise Exception('can not run docker command without exact version tag')
        image_name = self.image_name_with_tag.split(':', 1)[0]
        container_name = "{}_{}".format(image_name, datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))
        volumes = {self.navigomix_folder: {'bind': '/navigomix', 'mode': 'rw'},
                   self.input_folder: {'bind': '/input', 'mode': 'rw'},
                   self.output_folder: {'bind': '/output', 'mode': 'rw'}}

        if not self.silent:
            self.display.print("\n\n======= running docker image %s with command: %s (volumes: %s)" % (self.image_name_with_tag, command, str(volumes)))
        container = self.docker_client_high_level.containers.run(self.image_name_with_tag, command, working_dir='/', name=container_name,
                                                                 volumes=volumes, stdout=True, stderr=True, detach=True,
                                                                 environment={'PYTHONPATH': '/analytic-modules'})

        return container
