######## /usr/bin/python
# Ref: https://developers.google.com/youtube/1.0/developers_guide_python#UploadingVideos

import gdata.youtube
import gdata.youtube.service

yt_service = gdata.youtube.service.YouTubeService()

# Turn on HTTPS/SSL access.
# Note: SSL is not available at this time for uploads.
yt_service.ssl = True


####### Authentication 

# yt_service.client_id=''

yt_service.email='<EMAIL>@gmail.com'
yt_service.password='<PASSWORD>'
yt_service.source='ceor example'
yt_service.developer_key='AI39si6c1zCpV6asOC4hLdel_Df__enFJxHcgJrdnBrTcTyO81PLq0C7Q5KxVvGahdcGyCdBhyVlmrVqpdgytpZvV_b27g7TvA'
#yt_service.client_id='ceor_s example'
yt_service.ProgrammaticLogin()

my_media_group = gdata.media.Group(title=gdata.media.Title(text='My Test Movie'),description=gdata.media.Description(description_type='plain',text='My description'),keywords=gdata.media.Keywords(text='cars, funny'),category=gdata.media.Category(text='Autos',scheme='http://gdata.youtube.com/schemas/2007/categories.cat',label='Autos'),player=None )


where = gdata.geo.Where()
where.set_location((-22.979794,-43.232435))

# create the gdata.youtube.YouTubeVideoEntry to be uploaded
video_entry = gdata.youtube.YouTubeVideoEntry(media=my_media_group,
                                              geo=where)

# set the path for the video file binary
video_file_location = '/home/alfaceor/data/3D_1000100010001_0.04/tmp/untitled2.mpg'

new_entry = yt_service.InsertVideoEntry(video_entry, video_file_location)
