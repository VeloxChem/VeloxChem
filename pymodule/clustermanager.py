import os
import getpass
import pexpect
import sys
import time
import re

class ClusterManager:
    def __init__(self, cluster_name, path_on_cluster=None):
        self.cluster_name = cluster_name
        self.path_on_cluster_given = True
        self.job_id_file = f'{path_on_cluster}/slurm_job_id.txt'
        self.path_on_cluster = path_on_cluster
        if self.path_on_cluster is None:
            self.path_on_cluster_given = False
            self.path_on_cluster = os.getcwd()
            self.job_id_file = 'slurm_job_id.txt'
        self._password = None
        self.child = None
        
    def get_password(self):
        if self._password is None:
            self._password = getpass.getpass(f"Enter password for {self.cluster_name}: ")
            print("Password stored for reuse.")
        return self._password

    def connect(self, jobscript_path, input_files, geometry):
        if self.child is None or not self.child.isalive():
            if self.path_on_cluster_given is True:
                print('I am performing the scp command', jobscript_path)
                scp_command = f'scp {jobscript_path} {input_files} {geometry} {self.cluster_name}:{self.path_on_cluster}/.'
                self.child = pexpect.spawn(scp_command, encoding='utf-8', timeout=60)
                # self.child.logfile = sys.stdout  # Log all output to stdout for debugging
            
                i = self.child.expect(['[Pp]assword:', pexpect.TIMEOUT, pexpect.EOF], timeout=30)
                if i == 0:
                    self.child.sendline(self.get_password())
                i = self.child.expect(['[$#>]', pexpect.TIMEOUT, pexpect.EOF], timeout=30)
                print('I am performing the scp command', scp_command, i)
                self.child.sendline('exit')
                self.child.close()

            ssh_command = f'ssh -o ConnectTimeout=30 -o ServerAliveInterval=30 {self.cluster_name}'
            self.child = pexpect.spawn(ssh_command, encoding='utf-8', timeout=60)
            self.child.logfile = sys.stdout  # Log all output to stdout for debugging

            i = self.child.expect(['[Pp]assword:', pexpect.TIMEOUT, pexpect.EOF], timeout=30)
            if i == 0:
                self.child.sendline(self.get_password())
            elif i == 1:
                print("Timeout while waiting for password prompt.")
                return False
            elif i == 2:
                print("EOF encountered while waiting for password prompt.")
                return False

            i = self.child.expect(['\$ ', '# ', pexpect.TIMEOUT], timeout=60)
            if i == 0:
                print("Successfully connected to the cluster.")
                return True
            else:
                print("Failed to get shell prompt after password entry.")
                return False
        return True

    def disconnect(self):
        if self.child and self.child.isalive():
            self.child.sendline('exit')
            self.child.close()
            print('close child')
        self.child = None

    def run_cluster_command(self, command, jobscript_path, input_files, geometry):
        if not self.connect(jobscript_path, input_files, geometry):
            self.connect(jobscript_path, input_files, geometry)
            #return ""

        full_command = f"cd {self.path_on_cluster} && {command}"
        #print(f"Full command being executed: {full_command}")
        self.child.sendline(full_command)
        self.child.expect(['\n', pexpect.TIMEOUT, pexpect.EOF], timeout=5)  # Wait for command echo

        output = ""
        while True:
            i = self.child.expect(['\$', pexpect.EOF, pexpect.TIMEOUT], timeout=20)
            output += self.child.before  # Collect the text before the match
            if i == 0:  # Shell prompt detected
                break
            elif i == 1:  # EOF detected
                break
            elif i == 2:  # Timeout
                print("Timeout reached while waiting for command output.")
                break

        return output.strip()
    
    def extract_job(self, command):

        full_command = f"{command}"
        #print(f"Full command being executed: {full_command}")
        self.child.sendline("stty -echo")
        self.child.expect('\$ ')
        self.child.sendline(full_command)
        self.child.expect(['\$ ', '# ', pexpect.EOF], timeout=60)  # Wait for the prompt

        # Retrieve the output
        output = self.child.before.strip()  # Strip leading/trailing whitespace

        # Clean up known unwanted parts if needed
        if "Submitted batch job" in output:
            output = output.split("Submitted batch job")[-1].strip()  # Keep content after known phrase

        # Remove unexpected lines
        lines = output.split('\n')  # Split into lines
        filtered_lines = [line for line in lines if "rom_local" not in line and "LYRA" not in line]

        # Join the filtered lines back into a string
        cleaned_output = '\n'.join(filtered_lines)
    
        return cleaned_output

    def submit_job(self, jobscript_path, input_files, geometry, max_attempts=5, retry_delay=2):
        command = f'sbatch {jobscript_path} > {self.job_id_file}'
        output = self.run_cluster_command(command, jobscript_path, input_files, geometry)
        time.sleep(5)
        #print(f"Raw output from sbatch command: {repr(output)}")
        command = f"cat {self.job_id_file}"  # Command to read the file
        file_content = self.extract_job(command)
        print('FILE CONTENT', file_content)
        job_id = file_content
        print('JOB_ID', job_id)
        if job_id:
                print(f"Job submitted with ID: {job_id}")
                return job_id
        else:
            print(f"Failed to find job ID after {max_attempts} attempts. Check if the job was actually submitted.")
            exit()
            return None
        
        # attempts = 0

        # while attempts < max_attempts:
        #     job_id = self.read_job_id_from_file()
        #     if job_id:
        #         print(f"Job submitted with ID: {job_id}")
        #         return job_id
        #     else:
        #         print(f"Attempt {attempts + 1}/{max_attempts}: Couldn't find job ID in the output file. Retrying...")
        #         attempts += 1
        #         if attempts < max_attempts:
        #             time.sleep(retry_delay)

        # print(f"Failed to find job ID after {max_attempts} attempts. Check if the job was actually submitted.")
        # return None

    def read_job_id_from_file(self):
        try:
            # Read the file directly from the local filesystem
            files_in_current_directory = os.listdir('.')

            # Print the list of files
            for file_name in files_in_current_directory:
                print(file_name)
                
            with open(self.job_id_file, 'r') as file:
                content = file.read().strip()
            print(f"Content of {self.job_id_file}: {content}")

            # Extract the job ID (last number in the string)
            parts = content.split()
            if len(parts) >= 4 and parts[-4] == "Submitted" and parts[-3] == "batch" and parts[-2] == "job":
                job_id = parts[-1]
                if job_id.isdigit():
                    print(f"Extracted job ID: {job_id}")
                    return job_id

            print("No valid job ID found in the file content.")
            return None
        except FileNotFoundError:
            print(f"Error: The file {self.job_id_file} was not found.")
            return None
        except Exception as e:
            print(f"Error reading job ID file: {e}")
            return None

    def is_job_running(self, job_id):
        print(f'Checking status for job ID: {job_id}')
        command_2 = f'sacct -j {job_id} --parsable2'
        command = f'squeue --job {job_id} --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"'
        output = self.run_cluster_command(command, None, None, None)
        output_2 = self.run_cluster_command(command_2, None, None, None)
        # Split the output into lines
        lines = output_2.strip().split('\n')

        # Check if we have exactly two lines (header + job info)
        id_found = False
        for line in lines:
            print('######', line)
            if job_id in line:
                print('JOB ID found program is still running')
                id_found = True
            else:
                # Only header line present, job is not in the queue
                print(f"Job {job_id} is no longer in the queue (finished or not found).")

        return id_found

    def run_cluster_job_and_wait(self, jobscript_path, input=None, geometry=None, queue=False):
        
        job_id = self.submit_job(jobscript_path, input, geometry)

        if job_id and queue == False:
            while True:
                if not self.is_job_running(job_id):
                    print(f"Job {job_id} has finished running.")
                    break
                print(f"Job {job_id} is still running. Waiting...")
                time.sleep(5)  # Wait for 1 minute before checking again

        self.disconnect()
        time.sleep(10)
        return job_id
        print("Disconnected from the cluster.")


    def check_queue(self, job_ids):
        
        njobs = len(job_ids)
        while njobs > 0:
            for job_id in job_ids:
                if not self.is_job_running(job_id):
                    print(f"Job {job_id} has finished running.")
                    job_ids.remove(job_id)
                    njobs -= 1
            
            print(f"Job {job_ids} is still running. Waiting...")
            time.sleep(5)  # Wait for 1 minute before checking again

        self.disconnect()

        return njobs
        print("Disconnected from the cluster.")


    def connect_copy(self, outputfile):
        print('I am performing the scp command', outputfile)
        scp_command = f'scp {self.cluster_name}:{self.path_on_cluster}/{outputfile} .'
        self.child = pexpect.spawn(scp_command, encoding='utf-8', timeout=60)
        # self.child.logfile = sys.stdout  # Log all output to stdout for debugging
    
        i = self.child.expect(['[Pp]assword:', pexpect.TIMEOUT, pexpect.EOF], timeout=30)
        if i == 0:
            self.child.sendline(self.get_password())
            i = self.child.expect(['[$#>]', pexpect.TIMEOUT, pexpect.EOF], timeout=30)
        print('I am performing the scp command', scp_command, i)
        self.disconnect()

        