"""Manage Pages using Git's existing credential helper; never log credentials."""
import argparse
import json
import os
import subprocess
import urllib.error
import urllib.request


def request(repo, method='GET', body=None, endpoint=None):
    environment = dict(os.environ, GIT_TERMINAL_PROMPT='0', GCM_INTERACTIVE='never')
    credential = subprocess.run(['git', 'credential', 'fill'],
                                input='protocol=https\nhost=github.com\n\n',
                                capture_output=True, text=True, env=environment)
    fields = dict(line.split('=', 1) for line in credential.stdout.splitlines() if '=' in line)
    token = fields.get('password')
    if not token and method != 'GET':
        raise RuntimeError('Git credential helper returned no GitHub authentication')
    req = urllib.request.Request(endpoint or f'https://api.github.com/repos/treygreer/{repo}/pages',
                                 method=method, data=None if body is None else json.dumps(body).encode(),
                                 headers={'Accept': 'application/vnd.github+json',
                                          **({'Authorization': f'Bearer {token}'} if token else {}),
                                          'User-Agent': 'treb-pages-deployment',
                                          'X-GitHub-Api-Version': '2022-11-28'})
    try:
        with urllib.request.urlopen(req, timeout=30) as response:
            return json.load(response) if response.status != 204 else {}
    except urllib.error.HTTPError as error:
        # GitHub's API response does not include credentials. Report only its message.
        payload = json.loads(error.read())
        raise RuntimeError(f'GitHub API {error.code}: {payload.get("message", "request failed")}') from None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repo', choices=['treb', 'firstinfright'])
    parser.add_argument('--create', action='store_true')
    parser.add_argument('--branch')
    parser.add_argument('--account', action='store_true')
    parser.add_argument('--https', action='store_true', help='Enforce HTTPS for this Pages site')
    parser.add_argument('--build', action='store_true', help='Read the latest Pages build diagnostics')
    args = parser.parse_args()
    if args.build:
        build = request(args.repo, endpoint=f'https://api.github.com/repos/treygreer/{args.repo}/pages/builds/latest')
        print(json.dumps({key: build.get(key) for key in ['status', 'error', 'commit', 'created_at', 'updated_at']}, indent=2))
        runs = request(args.repo, endpoint=f'https://api.github.com/repos/treygreer/{args.repo}/actions/runs?per_page=3')
        print(json.dumps({'runs': [{key: run.get(key) for key in ['name', 'status', 'conclusion', 'html_url']}
                                  for run in runs.get('workflow_runs', [])]}, indent=2))
        active = runs.get('workflow_runs', [])
        if active and active[0].get('status') != 'completed':
            jobs = request(args.repo, endpoint=f'https://api.github.com/repos/treygreer/{args.repo}/actions/runs/{active[0]["id"]}/jobs')
            print(json.dumps({'jobs': [{'name': job.get('name'), 'status': job.get('status'),
                                      'steps': [{key: step.get(key) for key in ['name', 'status', 'conclusion']}
                                                for step in job.get('steps', [])]}
                                     for job in jobs.get('jobs', [])]}, indent=2))
        return
    if args.account:
        user = request(args.repo, endpoint='https://api.github.com/user')
        repo = request(args.repo, endpoint=f'https://api.github.com/repos/treygreer/{args.repo}')
        print(json.dumps({'login': user.get('login'), 'permissions': repo.get('permissions')}, indent=2))
        return
    if args.create:
        if args.repo != 'treb' or not args.branch:
            parser.error('Creation is supported only for the treb viewer with an explicit branch')
        try:
            request(args.repo, 'POST', {'source': {'branch': args.branch, 'path': '/'}})
        except RuntimeError as error:
            if not str(error).startswith('GitHub API 409:'):
                raise
            existing = request(args.repo)
            if existing.get('source') != {'branch': args.branch, 'path': '/'}:
                raise RuntimeError('Pages is already enabled with a different publishing source') from None
    if args.https:
        request(args.repo, 'PUT', {'https_enforced': True})
    data = request(args.repo)
    print(json.dumps({key: data.get(key) for key in
                      ['status', 'html_url', 'source', 'cname', 'https_enforced', 'build_type']}, indent=2))


if __name__ == '__main__':
    main()
