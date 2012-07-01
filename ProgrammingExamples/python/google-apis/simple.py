import gdata.docs.data
import gdata.docs.client

client = gdata.docs.client.DocsClient(source='myapp_ceor-v0')
client.ssl = True
client.http_client.debug = False
client.ClientLogin('emailr@gmail.com','password',client.source);

def PrintFeed(feed):
  print '\n'
  if not feed.entry:
    print 'No entries in feed.\n'
  for entry in feed.entry:
    print entry.title.text.encode('UTF-8'), entry.GetDocumentType(), entry.resource_id.text

    # List folders the document is in.
    for folder in entry.InFolders():
      print folder.title

feed = client.GetDocList(uri='/feeds/betaceor@gmail.com/private/full?showfolders=true')
PrintFeed(feed)

#new_folder= client.Create(gdata.docs.data.FOLDER_LABEL, 'Research')
import gdata.data
ms= gdata.data.MediaSource(file_path='/home/alfaceor/Research/Eskeleton_ProteinFolding.txt', content_type='text/plain')

#feed = client.GetDocList(uri='/feeds/default/private/full?title=Research&title-exact=true&max-results=5')
feed = client.GetDocList(uri='/feeds/default/private/full/-/folder')
PrintFeed(feed)
research_folder='https://docs.google.com/feeds/default/private/full/folder%3A0B_kALGdOkvLYOGY4MDA0NTItMzFiNS00YTA2LThiYzYtOTA2NjJlMTM5Njc1/contents'
entry=client.Upload(ms,'Eskeleton_ProteinFolding',folder_or_uri=research_folder)
