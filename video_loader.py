
from pathlib import Path
import pickle

class VideoLoader:
    def __init__(self, base_dir, metadata_file):
        self.base_dir = Path(base_dir)
        with open(self.base_dir / metadata_file, 'rb') as f:
            self.metadata = pickle.load(f)
    
    def get_video(self, cube, condition, line):
        try:
            video_path = self.base_dir / self.metadata[cube][condition][line]
            with open(video_path, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            print(f"Error loading video: {str(e)}")
            return None
    
    def list_videos(self):
        return [(cube, condition, line) 
                for cube in self.metadata
                for condition in self.metadata[cube] 
                for line in self.metadata[cube][condition]]