# /usr/bin/env python

import sys

def sing_command():
    return "singularity exec /gpfs/projects/bsc05/guidance/singularity_images/guidance_singularity_compss25_06_11_19" \
           ".img "


def java_jar():
    return "java -cp /gpfs/projects/bsc05/guidance/guidance_versions/guidance_25_01_02_21.jar " \
           "guidance.FunctionWrappers "


"""
def transform_command_to_singularity(command):
    command = command.replace("/usr/bin/plink",
                              sing_command() + "/usr/bin/plink")
    command = command.replace("/usr/lib/jvm/java-8-openjdk-amd64/java createRsIdList",
                              sing_command() + java_jar() + "3")
    command = command.replace("/TOOLS/shapeit.v2.r727.linux.x64",
                              sing_command() + "/TOOLS/shapeit.v2.r727.linux.x64")
    command = command.replace("/usr/lib/jvm/java-8-openjdk-amd64 newSample.jar", sing_command() + java_jar() + "7")
    command = command.replace("/usr/lib/jvm/java-8-openjdk-amd64/java imputeWithImputeAndFilterByInfo.jar",
                              sing_command() + java_jar() + "8")
    return command
"""

class Command:
    def __init__(self, line, counter):
        self.line = line
        self.counter = counter

    def set_commands(self, commands):
        self.commands = commands

    def full_command(self):
        return self.dependency() + " " + self.command()

    def dependency(self):
        return ""

    def line_for_string(self, string_list):
        for command in self.commands:
            """All the words passed in the list are present in the line"""
            if reduce(lambda x, y: x and y, map(lambda x: x in command.line, string_list)):
                return str(command.counter)
        return ""


class FromBedToBed(Command):

    def command(self):
        return self.line.replace("/usr/bin/plink",
                                 sing_command() + "/usr/bin/plink")


class CreateRsIdList(Command):

    def command(self):
        return self.line.replace("/usr/lib/jvm/java-8-openjdk-amd64/java createRsIdList",
                                 sing_command() + java_jar() + "3")

    def dependency(self):
        string_to_search = self.line.split()[2][:-4]
        return "[# " + self.line_for_string([string_to_search]) + " #]"


class Phasing(Command):

    def command(self):
        return self.line.replace("/TOOLS/shapeit.v2.r727.linux.x64",
                                 sing_command() + "/TOOLS/shapeit.v2.r727.linux.x64")

    def dependency(self):
        string_to_search = self.line.split()[2][:-4]
        return "[# " + self.line_for_string([string_to_search]) + " #]"


class NewSample(Command):

    def command(self):
        return self.line.replace("/usr/lib/jvm/java-8-openjdk-amd64 newSample.jar", sing_command() + java_jar() + "7")

    def dependency(self):
        string_to_search = self.line.split()[3]
        return "[# " + self.line_for_string([string_to_search]) + " #]"


class Imputation(Command):

    def command(self):
        return self.line.replace("/usr/lib/jvm/java-8-openjdk-amd64/java imputeWithImputeAndFilterByInfo.jar",
                                 sing_command() + java_jar() + "8")

    def dependency(self):
        return "[# " + self.line_for_string([self.line.split()[9], "createRsIdList"]) + "," + self.line_for_string(
            [self.line.split()[6], "newSample.jar"]) + " #]"


def get_instance_from_line(line, counter):
    if "/usr/bin/plink" in line and "make-bed" in line:
        return FromBedToBed(line, counter)
    if "createRsIdList" in line:
        return CreateRsIdList(line, counter)
    if "/TOOLS/shapeit.v2.r727.linux.x64" in line:
        return Phasing(line, counter)
    if "newSample.jar" in line:
        return NewSample(line, counter)
    if "imputeWithImputeAndFilterByInfo.jar" in line:
        return Imputation(line, counter)
    raise ("ERROR")


def translations():
    return [
        ("/usr/bin/plink", sing_command() + "/usr/bin/plink"),
        ("/usr/lib/jvm/java-8-openjdk-amd64/java createRsIdList",
         sing_command() + java_jar() + "3"),
        ("/TOOLS/shapeit.v2.r727.linux.x64",
         sing_command() + "/TOOLS/shapeit.v2.r727.linux.x64"),
        ("/usr/lib/jvm/java-8-openjdk-amd64 newSample.jar", sing_command() + java_jar() + "7"),
        ("/usr/lib/jvm/java-8-openjdk-amd64/java imputeWithImputeAndFilterByInfo.jar",
         sing_command() + java_jar() + "8")
    ]


def line_contains_translation(line):
    for pair in translations():
        if pair[0] in line:
            return True
    return False


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    counter = 1
    lines = []
    with open(input_file) as infile:
        for line in infile:
            if not line[0] == "#" and line_contains_translation(line) and line.strip():
                lines.append(get_instance_from_line(line, counter))
                counter += 1
    with open(output_file, "w") as outfile:
        for line in lines:
            line.set_commands(lines)
            outfile.write(line.full_command())
