#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utility for downloading Twitter feed to CSV and JSON files.

By Jason Miller
"""

import sys
import json
import csv
import twitter
from NotForGitHub import ACCESS_TOKEN_KEY, ACCESS_TOKEN_SECRET, CONSUMER_KEY, CONSUMER_SECRET

def get_tweets(api=None, hash_tag=None, goal=100):
    """Access Twitter API using python-twitter library.

    The python-twitter library is here:
    https://github.com/bear/python-twitter
    https://python-twitter.readthedocs.io

    Uses this library call:
    GetSearch(self, term=None, raw_query=None,
        geocode=None, since_id=None, max_id=None,
        until=None, since=None, count=15, lang=None,
        locale=None, result_type='mixed',
        include_entities=None, return_json=False)
        
    Following examples in the library sample code,
    this asks for a few tweets at a time from recent to old.
    """

    bucket_size = 10   # debugging
    bucket_size = 100   # production
    all_tweets = api.GetSearch(count=bucket_size,term=hash_tag,lang="en")
    if all_tweets and isinstance(all_tweets,list) and len(all_tweets)>0:
        num_total = len(all_tweets)
        while num_total < goal:
            earliest_id = min(all_tweets, key=lambda x: x.id).id
            earliest_id = earliest_id - 1
            print("{} so far. Get tweets before {}".format(len(all_tweets),earliest_id))
            more_tweets = api.GetSearch(count=bucket_size,term=hash_tag,lang="en",max_id=earliest_id)
            if not more_tweets or len(more_tweets)==0:
                break
            all_tweets += more_tweets
            num_total = len(all_tweets)
    else:
        all_tweets=[]
    return all_tweets

def write_screen(timeline):
    """Write a list of tweets to screen for debugging."""
    for tweet in timeline:
        print ("Lan={} Len={}\n{}\n".format(tweet.lang,
            len(tweet.full_text),tweet.full_text.replace('\n','')))

def write_json(timeline,filename):
    """Write a list of tweets to a JSON file for browsing."""
    with open(filename, 'w') as f:
        for tweet in timeline:
            f.write(json.dumps(tweet._json,indent=4))
            f.write('\n')

def safe_get(dictionary,field,default):
    if not isinstance(dictionary,dict) or field is None or not field in dictionary:
        return default
    return dictionary.get(field)

def write_csv(timeline,filename,includeRetweet=False):
    """Write a list of tweets to a CSV file for processing."""
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        fieldnames = ['id_str','created_at','screen_name',
            'retweet_count','favorite_count','full_text']
        writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_NONNUMERIC)
        writer.writeheader()
        for tweet in timeline:
            one_tweet=tweet.AsDict()
            t_id = safe_get(one_tweet,'id_str','None')
            t_created = safe_get(one_tweet,'created_at','None')
            t_user = safe_get(one_tweet,'user','None')
            t_name = safe_get(t_user,'screen_name','None')
            t_retweet = safe_get(one_tweet,'retweet_count',0)
            t_favorite = safe_get(one_tweet,'favorite_count',0)
            t_text = safe_get(one_tweet,'full_text','None').replace('\n',' ')
            if includeRetweet or not t_text.startswith("RT "):
                sub_dict={'id_str':t_id,'created_at':t_created,'screen_name':t_name,
                    'retweet_count':t_retweet,'favorite_count':t_favorite,
                    'full_text':t_text}
                writer.writerow(sub_dict)

def twitter_initialize():
    """Prepare to use the python-twitter library."""
    api = twitter.Api(
        CONSUMER_KEY, CONSUMER_SECRET, ACCESS_TOKEN_KEY, ACCESS_TOKEN_SECRET,
        tweet_mode='extended'
    )
    return api
